(ns bcbio.variation.ensemble.realign
  "Realignment based Ensemble calling approaches using inputs from multiple callers.
   Uses tools from the Marth lab and Erik Garrison to realign and recall given a
   set of possible variants in a genomic region."
  (:require [bcbio.variation.ensemble.prep :as eprep]))

(defn by-region
  "Realign and recall variants in a defined genomic region."
  [region vcf-files bam-file ref-file]
  (let [prep-inputs (map eprep/norm-bgzip vcf-files)]
    (println prep-inputs)))
