(ns bcbio.variation.ensemble.prep
  "Prepare variant inputs for ensemble calling approaches.
   Creates normalized, bgzipped tabix indexed inputs from VCF files."
  (:require [bcbio.run.itx :as itx]
            [clojure.core.strint :refer [<<]]
            [me.raynes.conch.low-level :as conch]))

(defn check-run
  [cmd]
  (let [proc (conch/proc "bash" "-c" (str "set -o pipefail; " cmd))]
    (future (conch/stream-to-out proc :out))
    (future (conch/stream-to-out proc :err))
    (when-not (= 0 (conch/exit-code proc))
      (throw (Exception. (format "Shell command failed: %s" cmd))))))

(defn norm-bgzip
  "Normalize and bgzip/tabix index a VCF input file."
  [vcf-file]
  (let [out-file (str vcf-file ".gz")
        tabix-file (str out-file ".tbi")]
    (when (itx/needs-run? out-file)
      (itx/with-tx-file [tx-out-file out-file]
        (check-run (<< "vcfallelicprimitives ~{vcf-file} | "
                       "bgzip -c /dev/stdin > ~{tx-out-file}"))))
    (when (itx/needs-run? tabix-file)
      (check-run (<< "tabix -p vcf ~{out-file}")))
    out-file))
