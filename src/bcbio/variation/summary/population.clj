(ns bcbio.variation.summary.population
  "Summarize variant metrics and extract useful variants for batch called populations."
  (:import [net.sf.samtools SAMFileReader]
           [net.sf.picard.util SamLocusIterator IntervalList Interval])
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [criterium.stats :as stats]
            [doric.core :as doric]
            [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Filtering

(defn- passes-filter?
  "Does variant pass the original GATK filter."
  [vc]
  (empty? (:filters vc)))

(defn- passes-call-rate?
  "Check for 100% variant call rate."
  [vc]
  (not (contains? (set (map :type (:genotypes vc))) "NO_CALL")))

(defn- has-min-depth?
  "Ensure all samples have at least 12x coverage for calling heterozygotes."
  [depth vc]
  (every? #(>= % depth) (map #(get-in % [:attributes "DP"]) (:genotypes vc))))

(defn- has-variable-genotypes?
  "Check that a genotype varies across the experimental conditions."
  [vc]
  (> (count (set (map :type (:genotypes vc))))
     1))

(defn- het-in-one-treatment?
  "Ensure variant is heterozygous in only one treatment condition"
  [samples vc]
  (let [treats (reduce (fn [coll g]
                         (if (= (:type g) "HET")
                           (conj coll (get samples (:sample-name g)))
                           coll))
                       #{} (:genotypes vc))]
    (= 1 (count treats))))

(defn- get-target-depth
  "Retrieve desired depth for a set of samples from configuration file.
   Defaults to heterozygous resolution depth (12)"
  [vcf-file config-file]
  (let [config (-> config-file slurp yaml/parse-string)
        sample-names (-> vcf-file gvc/get-vcf-header .getGenotypeSamples vec)
        depths (reduce (fn [coll x]
                         (if-let [depth (get-in x [:metadata :depth])]
                           (assoc coll (:description x) depth)
                           coll))
                       {} (:details config))
        cur-depths (vals (select-keys depths sample-names))]
    (if (empty? cur-depths) 12 (apply min cur-depths))))

(defn- get-bam-depth
  "Attempt at lightweight way to get depth at a position."
  [bam-file chrom start]
  (with-open [rdr (SAMFileReader. (io/file bam-file) (io/file (str bam-file ".bai")))
              sli (.iterator (SamLocusIterator. rdr (doto (IntervalList. (.getFileHeader rdr))
                                                      (.add (Interval. chrom (inc start) (inc start))))
                                                true))]
    (-> sli .next .getRecordAndPositions .size)))

(defn check-evaluation-depth
  "Prepare function that ensures a set of baseline BAM files have minimum depth at a position"
  [bam-files]
  (let [depth 4]
    (fn [vc]
      (every? #(>= (get-bam-depth % (:chr vc) (:start vc)) depth) bam-files))))

;; ## Localize variants per gene

(defn- parse-snpeff-line
  [vc]
  (letfn [(eff->gene-name [eff]
            (let [priorities {"HIGH" 0 "MODERATE" 1 "LOW" 2 "MODIFIER" 3}
                  parts (-> eff
                            (string/split #"\(")
                            second
                            (string/split #"\|"))
                  size (nth parts 4)]
              {:gene (nth parts 5)
               :change (first parts)
               :priority (get priorities (first parts) 4)
               :codon-change (nth parts 2)
               :aa-change (nth parts 3)
               :size (when-not (empty? size) (* 3 (Integer/parseInt size)))}))]
    (map eff->gene-name
       (let [a (get-in vc [:attributes "EFF"])]
         (if (string? a) [a] a)))))

(defn- vc->gene-names
  "Extract gene names and sizes from snpEff annotation of variant effects.."
  [vc]
  (reduce (fn [coll g]
            (if (:size g)
              (assoc coll (:gene g) (:size g))
              coll))
          {} (parse-snpeff-line vc)))

(defn- update-gene-count
  "Extract gene names effected by variant and update count."
  [orig-coll vc]
  (reduce (fn [coll [g size]]
            (assoc coll g
                   (update-in (get coll g {:count 0 :size size}) [:count] inc)))
          orig-coll (vc->gene-names vc)))

(defn- get-top-genes
  "Retrieve the top genes with large numbers of mutations, looking for outliers."
  [gcounts]
  (let [counts (sort > (map second gcounts))
        thresh (stats/quantile 0.3 counts)]
    (->> gcounts
         (sort-by second >)
         (take-while #(> (second %) thresh))
         (into {}))))

(defn- normalize-counts
  "Normalize variant counts by gene size"
  [gcounts]
  (reduce (fn [coll [g c]]
            (assoc coll g (/ (:count c) (:size c))))
          {} gcounts))

(defn- low-count?
  [[_ g]]
  (< (:count g) 5))

(defn highly-mutated-genes
  "Identify genes with large numbers of mutations for filtering purposes."
  [vcf-file config-file]
  (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file)]
    (->> (gvc/parse-vcf vcf-iter)
         (filter (partial has-min-depth? (get-target-depth vcf-file config-file)))
         (reduce update-gene-count {})
         (remove low-count?)
         normalize-counts
         get-top-genes)))

(defn check-highly-mutated
  "Generate predicate function to identify variants in highly mutated genes."
  [vcf-file input-genes config-file]
  (let [high-genes (highly-mutated-genes vcf-file config-file)
        check-genes (set/union (set (map #(string/lower-case (first %)) high-genes)) input-genes)]
    (println "** Highly mutated")
    (doseq [[k v] high-genes]
       (println k v))
    (fn [vc]
      (some (fn [[k v]] (contains? check-genes (string/lower-case k)))
            (vc->gene-names vc)))))

(defn get-clustered-genes
  "Retrieve a list of manually identified clustered genes to avoid in summarizing."
  [fname]
  (letfn [(find-gene-coords [parts]
            (->> parts
                 (map-indexed (fn [i x] (when (= x "gene") i)))
                 (remove nil?)))
          (get-genes [parts coords]
            (when (seq (filter #(.startsWith % "chr") parts))
              (->> coords
                   (map #(nth parts %))
                   (map string/lower-case)
                   (remove empty?)
                   set)))]
    (with-open [reader (io/reader fname)]
      (-> (reduce (fn [coll parts]
                    (cond
                     (= "chr" (first parts)) (assoc coll :gene-coords (find-gene-coords parts))
                     (:gene-coords coll) (assoc coll :genes (set/union (:genes coll)
                                                                       (get-genes parts (:gene-coords coll))))
                     :else coll))
                  {:gene-coords nil :genes #{}} (csv/read-csv reader))
          :genes))))

;; ## Organize samples and treatments

(defn- get-sample-treatments
  "Retrieve map of sample names to treatments from an input YAML metadata."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (reduce (fn [coll x]
              (assoc coll (:description x) (get-in x [:metadata :treatment])))
            {} (:details config))))

(defn- get-treatment-map
  "Retrieve map of sample names to treatments."
  [vcf-file config-file]
  (let [sample-names (-> vcf-file gvc/get-vcf-header .getGenotypeSamples vec)]
    (select-keys (get-sample-treatments config-file) sample-names)))

;; ## Summarize calls by treatment

(defn- vc->treat-types
  "Retrieve the treatment and variant call types for a variant"
  [samples vc]
  (let [types-by-treat (reduce (fn [coll g]
                                 (let [k (get samples (:sample-name g))]
                                   (assoc coll k (conj (get coll k #{}) (:type g)))))
                               {} (:genotypes vc))]
    (map (fn [[treat types]]
              (let [val (cond
                         (contains? types "HET") "HET"
                         (contains? types "HOM_REF") "HOM_REF"
                         (contains? types "HOM_VAR") "HOM_VAR")]
                [treat (keyword val)]))
         types-by-treat)))

(defn- vc-treatment-calls
  "Summarize calls by treatment for a variant context."
  [samples orig-coll vc]
  (reduce (fn [coll key]
            (assoc-in coll key (inc (get-in coll key 0))))
          orig-coll (vc->treat-types samples vc)))

(defn- print-metrics
  [samples calls-by-treatment]
  (letfn [(summarize-treat [calls-by-treatment treat]
            (let [total (apply + (vals (get calls-by-treatment treat)))
                  treat-samples (reduce (fn [coll [k v]]
                                          (assoc coll v (conj (get coll v []) k)))
                                        {} samples)
                  treat-str (format "%s (%s)" treat (count (get treat-samples treat)))]
              (reduce (fn [coll x]
                        (let [count (get-in calls-by-treatment [treat x] 0)]
                          (-> coll
                              (assoc x count)
                              (assoc (keyword (str (name x) "-%"))
                                     (format "%.1f" (* 100.0 (/ count total)))))))
                      {:treatment treat-str} [:HOM_REF :HET :HOM_VAR])))]
    (println (doric/table [:treatment :HOM_REF :HOM_REF-% :HET :HET-% :HOM_VAR :HOM_VAR-%]
                          (map (partial summarize-treat calls-by-treatment) (sort (keys calls-by-treatment)))))))

(defn- print-chromosome-distribution
  "Table of reads mapped by chromosome."
  [vcs]
  (letfn [(gs->het-counts [top-coll vc]
            (reduce (fn [coll g]
                      (if (= "HET" (:type g))
                        (assoc coll (:sample-name g)
                               (inc (get coll (:sample-name g) 0)))
                        coll))
                    top-coll (:genotypes vc)))
          (chr->sample-counts [chrom samples counts]
            (reduce (fn [coll s]
                      (assoc coll s (get-in counts [chrom s] 0)))
                    {:chr chrom} samples))]
    (let [counts (reduce (fn [coll vc]
                           (assoc coll (:chr vc)
                                  (gs->het-counts (get coll (:chr vc)) vc)))
                         {} vcs)
          samples (map :sample-name (:genotypes (first vcs)))]
      (println (doric/table (cons :chr samples)
                            (map #(chr->sample-counts % samples counts) (sort (keys counts)))))))
  vcs)

(defn call-metrics
  "Calculate overall summary call metrics for an input VCF file."
  [vcf-file config-file checks]
  (println "**" (fs/base-name vcf-file))
  (let [samples (get-treatment-map vcf-file config-file)]
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file)]
      (->> (gvc/parse-vcf vcf-iter)
           (filter passes-filter?)
           (filter passes-call-rate?)
           (filter (partial has-min-depth? (get-target-depth vcf-file config-file)))
           (filter has-variable-genotypes?)
           ;(filter (partial het-in-one-treatment? samples))
           (remove (:is-high-gene? checks))
           (filter (:passes-depth? checks))
           print-chromosome-distribution
           (reduce (fn [coll vc] (vc-treatment-calls samples coll vc)) {})
           (print-metrics samples)))))

;; ## Summary CSV of calls

(defn- snpeff-effect
  "Summarize most worrisome effect from snpEff output."
  [vc]
  (let [e (first (sort-by :priority (parse-snpeff-line vc)))]
    (map #(get e %) [:gene :change :codon-change :aa-change])))

(defn- vc->summary
  "Covert a variant into a single line CSV-friendly summary."
  [samples vc]
  (concat
   [(:chr vc)
    (:start vc)
    (.getBaseString (:ref-allele vc))
    (string/join "&" (map #(.getBaseString %) (:alt-alleles vc)))
    (int (stats/mean (map #(get-in % [:attributes "DP"]) (:genotypes vc))))
    (string/join "&" (map first (filter #(= (second %) :HET) (vc->treat-types samples vc))))]
   (snpeff-effect vc)
   (map #(:type %) (:genotypes vc))))

(defn call-summary-csv
  "Summarize calls into a final CSV for exploration."
  [vcf-file config-file checks]
  (let [samples (get-treatment-map vcf-file config-file)
        out-file (str (fsp/file-root vcf-file) "-summary.csv")]
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file)
                wtr (io/writer out-file)]
      (csv/write-csv wtr [(concat ["chr" "start" "ref" "alt" "depth" "treatment" "gene" "effect" "codon" "aa"]
                                  (-> vcf-file gvc/get-vcf-header .getGenotypeSamples))])
      (->> (gvc/parse-vcf vcf-iter)
           (filter passes-call-rate?)
           (filter (partial has-min-depth? (get-target-depth vcf-file config-file)))
           (filter has-variable-genotypes?)
           ;(filter (partial het-in-one-treatment? samples))
           (remove (:is-high-gene? checks))
           (filter (:passes-depth? checks))
           (map (partial vc->summary samples))
           (csv/write-csv wtr)))))

;; ## Evaluate replicate consistency

(defn- evalute-consistency
  [samples orig-coll vc]
  (let [calls-by-treat
        (reduce (fn [calls g]
                  (let [treat (get samples (:sample-name g))]
                    (assoc calls treat (conj (get calls treat #{}) (:type g)))))
                {} (:genotypes vc))]
    (reduce (fn [coll [treat calls]]
              (let [kw (if (= 1 (count calls)) :consistent :multiple)]
                (assoc-in coll [treat kw] (inc (get-in coll [treat kw] 0)))))
            orig-coll calls-by-treat)))

(defn- print-consistency
  [samples summary]
  (let [treat-samples (reduce (fn [coll [k v]]
                                (assoc coll v (conj (get coll v []) k)))
                              {} samples)]
    (doseq [treat (sort (keys summary))]
      (println "***" treat (count (get treat-samples treat)))
      (let [total (apply + (vals (get summary treat)))]
        (doseq [con-type (sort (keys (get summary treat)))]
          (let [cur-count (get-in summary [treat con-type])]
            (println (format " - %s %s %.1f%%" (name con-type) cur-count
                             (* 100.0 (/ cur-count total))))))))))

(defn consistency-metrics
  "Calculate metrics of consistency of calls amongst replicates."
  [vcf-file config-file checks]
  (let [samples (get-treatment-map vcf-file config-file)]
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file)]
      (->> (gvc/parse-vcf vcf-iter)
           (filter passes-call-rate?)
           (filter (partial has-min-depth? (get-target-depth vcf-file config-file)))
           (filter has-variable-genotypes?)
           (remove (:is-high-gene? checks))
           (filter (:passes-depth? checks))
           (reduce (partial evalute-consistency samples) {})
           (print-consistency samples)))))

(defn- input-files-by-ext
  "Split input files by file extensions."
  [fs]
  (reduce (fn [coll f]
            (let [ext (fs/extension f)]
              (assoc coll ext
                     (cons f (get coll ext [])))))
          {} fs))

(defn call-metrics-many
  [config-file cluster-file & support-files]
  (let [{bam-files ".bam" vcf-files ".vcf"} (input-files-by-ext support-files)
        clustered-genes (get-clustered-genes cluster-file)]
    (doseq [vcf-file vcf-files]
      (let [checks {:is-high-gene? (check-highly-mutated vcf-file clustered-genes config-file)
                    :passes-depth? (check-evaluation-depth bam-files)}]
        (call-metrics vcf-file config-file checks)
        (call-summary-csv vcf-file config-file checks)
        (consistency-metrics vcf-file config-file checks)))))

(defn -main
  [& args]
  (if (>= (count args) 3)
    (apply call-metrics-many args)
    (do
      (println "Usage:")
      (println "  pop-summary <config-file> <cluster-file> <vcf-and-bam-files>"))))

(defn- tester
  []
  (let [;work-dir "/home/chapmanb/tmp/vcf/mouse"
        work-dir "/home/bchapman/tmp/dr_mouseexome/vcf"
        config-file (str work-dir "/drme-pilot.yaml")
        cluster-file (str work-dir "/clustered-paula.csv")
        vcf-file (str work-dir "/old1-gatk.vcf")
        bam-file (str work-dir "/1_130622_drme-sort.bam")]
    (call-metrics-many config-file cluster-file vcf-file bam-file)))
