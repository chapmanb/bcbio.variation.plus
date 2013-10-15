(ns bcbio.variation.summary.population
  "Summarize variant metrics and extract useful variants for batch called populations."
  (:require [clj-yaml.core :as yaml]
            [doric.core :as doric]
            [me.raynes.fs :as fs]
            [bcbio.variation.variantcontext :as gvc]))

(defn- get-sample-treatments
  "Retrieve map of sample names to treatments from an input YAML metadata."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (reduce (fn [coll x]
              (assoc coll (:description x) (get-in x [:metadata :treatment])))
            {} (:details config))))

(defn- passes-call-rate?
  "Check for 100% variant call rate."
  [vc]
  (not (contains? (set (map :type (:genotypes vc))) "NO_CALL")))

(defn- has-min-depth?
  "Ensure all samples have at least 12x coverage for calling heterozygotes."
  [vc]
  (every? #(>= % 12) (map #(get-in % [:attributes "DP"]) (:genotypes vc))))

(defn- vc-treatment-calls
  "Summarize calls by treatment for a variant context."
  [samples orig-coll vc]
  (reduce (fn [coll g]
            (let [key [(get samples (:sample-name g)) (keyword (:type g))]]
              (assoc-in coll key (inc (get-in coll key 0)))))
          orig-coll (:genotypes vc)))

(defn- print-metrics
  [calls-by-treatment]
  (letfn [(summarize-treat [calls-by-treatment treat]
            (let [total (apply + (vals (get calls-by-treatment treat)))]
              (reduce (fn [coll x]
                        (let [count (get-in calls-by-treatment [treat x] 0)]
                          (-> coll
                              (assoc x count)
                              (assoc (keyword (str (name x) "-%"))
                                     (format "%.1f" (* 100.0 (/ count total)))))))
                      {:treatment treat} [:HOM_REF :HET :HOM_VAR])))]
    (println (doric/table [:treatment :HOM_REF :HOM_REF-% :HET :HET-% :HOM_VAR :HOM_VAR-%]
                          (map (partial summarize-treat calls-by-treatment) (sort (keys calls-by-treatment)))))))

(defn call-metrics
  "Calculate overall summary call metrics for an input VCF file."
  [vcf-file ref-file config-file]
  (println "**" (fs/base-name vcf-file))
  (let [sample-names (-> vcf-file gvc/get-vcf-header .getGenotypeSamples vec)
        samples (select-keys (get-sample-treatments config-file) sample-names)]
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file ref-file)]
      (->> (gvc/parse-vcf vcf-iter)
           (filter has-min-depth?)
           (reduce (fn [coll vc] (vc-treatment-calls samples coll vc)) {})
           print-metrics))))

(defn- call-metrics-many
  [config-file ref-file & vcf-files]
  (doseq [vcf-file vcf-files]
    (call-metrics vcf-file ref-file config-file)))

(defn -main
  [& args]
  (if (>= (count args) 3)
    (apply call-metrics-many args)
    (do
      (println "Usage:")
      (println "  pop-summary <config-file> <ref-file> <vcf-files>"))))

(defn- tester
  []
  (let [ref-file "/usr/local/share/bcbio_nextgen/genomes/mm10/seq/mm10.fa"
        config-file "/home/chapmanb/tmp/vcf/mouse/drme-pilot.yaml"
        vcf-file "/home/chapmanb/tmp/vcf/mouse/old1-gatk.vcf"]
    (call-metrics vcf-file ref-file config-file)))
