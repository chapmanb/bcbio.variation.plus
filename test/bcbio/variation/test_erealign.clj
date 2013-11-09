(ns bcbio.variation.erealign
  "Tests for ensemble variant consolidation with local realignment."
  (:require [bcbio.variation.ensemble.realign :as erealign]
            [midje.sweet :refer :all]
            [clojure.java.io :as io]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (io/file "." "test" "data" "ensemble"))
               ref-file (str (io/file data-dir "chr10-start.fa"))
               bam-file (str (io/file data-dir "NA12878-10.bam"))
               vcf-files (map #(str (io/file data-dir %))
                              ["NA12878-10-freebayes.vcf" "NA12878-10-gatk.vcf"
                               "NA12878-10-gatk-haplotype.vcf"])
               work-dir (str (io/file data-dir "work"))]
           (doseq [x (concat [work-dir]
                             (map #(str % ".gz") vcf-files)
                             (map #(str % ".gz.tbi") vcf-files))]
             (itx/remove-path x))
           ?form)))

(facts "Calculate ensemble set of variants from multiple inputs using realignment."
  (let [region {:chrom "10" :start 250000 :end 400000}]
    (erealign/by-region region vcf-files bam-file ref-file work-dir)))
