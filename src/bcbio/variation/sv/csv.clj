(ns bcbio.variation.sv.csv
  "Extract structural variations from VCF into flat CSV structure"
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [bcbio.run.itx :as itx]
            [bcbio.variation.structural :as structural]
            [bcbio.variation.variantcontext :as gvc]))

(defn- vc->sv
  "Convert structural variation contexts into simple map structure"
  [vc]
  (when-let [sv-type (structural/get-sv-type vc {})]
    {:type sv-type
     :chr (:chr vc)
     :start (:start vc)
     :end (+ (:start vc) (#'structural/get-sv-length (assoc vc :sv-type sv-type)))
     :copy-number (get-in vc [:attributes "CN"] 1)
     :samples (->> (:genotypes vc)
                   (remove #(= "NO_CALL" (:type %)))
                   (map :sample-name))}))

(defn- sv->csv
  "Convert simple map structural variant into flat CSV output"
  [sv]
  (map (fn [sample]
         [sample (:chr sv) (:start sv) (:end sv) (name (:type sv)) (:copy-number sv)])
       (:samples sv)))

(defn vcf-convert
  "Convert structural variants in the input file to a flattened CSV."
  [vcf-file ref-file]
  (let [out-file (str (itx/file-root vcf-file) ".csv")]
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file ref-file)
                wtr (io/writer out-file)]
      (csv/write-csv wtr [["sample" "chr" "start" "end" "cnvtype" "cn"]])
      (csv/write-csv wtr (->> (gvc/parse-vcf vcf-iter)
                              (map vc->sv)
                              (remove nil?)
                              (mapcat sv->csv))))
    out-file))

(defn -main
  [& args]
  (if (= (count args) 2)
    (apply vcf-convert args)
    (do
      (println "Usage:")
      (println "  sv-to-csv <vcf-file> <ref-file>"))))
