(ns bcbio.variation.ensemble.realign
  "Realignment based Ensemble calling approaches using inputs from multiple callers.
   Uses tools from the Marth lab and Erik Garrison to realign and recall given a
   set of possible variants in a genomic region."
  (:require [bcbio.variation.ensemble.prep :as eprep]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn by-region
  "Realign and recall variants in a defined genomic region."
  [region vcf-files bam-file ref-file work-dir]
  (let [region-dir (io/file work-dir (eprep/region->safestr region))
        _ (when-not (fs/exists? region-dir) (fs/mkdirs region-dir))
        prep-inputs (map #(eprep/norm-bgzip % region region-dir) vcf-files)
        union-vcf (eprep/create-union prep-inputs ref-file region region-dir)]
    (println union-vcf)))
