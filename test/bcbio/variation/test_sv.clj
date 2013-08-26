(ns bcbio.variation.test-sv
  (:require [clojure.java.io :as io]
            [midje.sweet :refer :all]
            [bcbio.run.itx :as itx]
            [bcbio.variation.sv.csv :as sv.csv]))
(background
 (around :facts
         (let [data-dir (str (io/file "." "test" "data"))
               ref-file (str (io/file data-dir "GRCh37.fa"))
               sv-i-vcf (str (io/file data-dir "sv-illumina.vcf"))
               sv-i-vcf-out (str (itx/file-root sv-i-vcf) ".csv")]
           (doseq [x [sv-i-vcf-out]]
             (itx/remove-path x))
           ?form)))

(facts "Convert structural variants into flattened CSV structure"
  (sv.csv/vcf-convert sv-i-vcf ref-file) => sv-i-vcf-out)
