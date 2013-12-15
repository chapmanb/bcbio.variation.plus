(ns bcbio.variation.customrun.drme
  "Custom code for running top level analysis on mouse exome project."
  (:require [bcbio.variation.summary.population :as pop]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn run-summary
  []
  (let [ref-file "/n/hsphS10/hsphfs1/chb/biodata/genomes/Mmusculus/mm10/seq/mm10.fa"
        work-dir "/n/hsphS10/hsphfs1/chb/projects/dr_mouse_exomeseq_work/drme/"
        config-file (str work-dir "config/drme-pilot.yaml")
        vcf-files (sort (map str (fs/glob (str work-dir "exome_only/*.vcf"))))]
    (binding [*out* (io/writer "drme-summary.txt")]
      (apply pop/call-metrics-many (concat [config-file ref-file] vcf-files)))))
