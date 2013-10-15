(ns bcbio.variation.plus
  "Main namespace for exposing command line based functionality in bcbio.variation.plus version"
  (:require [bcbio.variation.summary.population]
            [bcbio.variation.sv.csv])
  (:gen-class))

(def ^{:private true} progs
  {:sv-to-csv 'bcbio.variation.sv.csv
   :pop-summary 'bcbio.variation.summary.population})

(defn -main [& args]
  (if-let [to-run-ns (get progs (keyword (first args)))]
    (apply (ns-resolve to-run-ns (symbol "-main"))
           (rest args))
    (do
      (println "Available commands:")
      (doseq [k (sort (keys progs))]
        (println (format " %s   %s" (name k) (:doc (meta (find-ns (get progs k))))))))))
