(ns bcbio.variation.ensemble.prep
  "Prepare variant inputs for ensemble calling approaches.
   Creates normalized, bgzipped tabix indexed inputs from VCF files."
  (:require [bcbio.run.itx :as itx]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.conch.low-level :as conch]
            [me.raynes.fs :as fs]))

(defn check-run
  [cmd]
  (let [proc (conch/proc "bash" "-c" (str "set -o pipefail; " cmd))]
    (future (conch/stream-to-out proc :out))
    (future (conch/stream-to-out proc :err))
    (when-not (= 0 (conch/exit-code proc))
      (throw (Exception. (format "Shell command failed: %s" cmd))))))

(defmacro run-cmd
  "Run a command line producing the given output file in an idempotent transaction.
   Wraps all of the machinery around preparing a command line from local arguments.
   Ensures a single run and avoids partial output files."
  [out-file & cmd]
  `(when (itx/needs-run? ~out-file)
     (itx/with-tx-file [tx-out-file# ~out-file]
       (let [fill-cmd# (<< ~@cmd)
             tx-cmd# (string/replace fill-cmd# ~out-file tx-out-file#)]
         (check-run tx-cmd#)))))

(defn- tabix-index-vcf
  [bgzip-file]
  (let [tabix-file (str bgzip-file ".tbi")]
    (run-cmd tabix-file "tabix -p vcf ~{bgzip-file}")
    tabix-file))

(defn bgzip-index-vcf
  "Prepare a VCF file for positional query with bgzip and tabix indexing."
  [vcf-file]
  (let [out-file (str vcf-file ".gz")]
    (run-cmd out-file "bgzip -c ~{vcf-file} > ~{out-file}")
    (tabix-index-vcf out-file)
    out-file))

(defn region->samstr
  [region]
  (format "%s:%s-%s" (:chrom region) (:start region) (:end region)))

(defn region->safestr
  [region]
  (format "%s_%s_%s" (:chrom region) (:start region) (:end region)))

(defn norm-bgzip
  "Normalize and bgzip/tabix index a VCF input file in a defined region."
  [vcf-file region out-dir]
  (let [prep-vcf-file (bgzip-index-vcf vcf-file)
        out-file (str (io/file out-dir (str (fs/base-name vcf-file) ".gz")))]
    (run-cmd out-file
             "tabix -h ~{prep-vcf-file} ~{(region->samstr region)} | "
             "vcfallelicprimitives | "
             "bgzip -c /dev/stdin > ~{out-file}")
    (tabix-index-vcf out-file)
    out-file))

(defn create-union
  "Create a union file with inputs from multiple variant callers."
  [vcf-files ref-file region out-dir]
  (let [out-file (str (io/file out-dir (str "union-" (region->safestr region) ".vcf")))
        intersect-str (string/join " | " (map (fn [x] (<< "vcfintersect -r ~{ref-file} -u ~{x}"))
                                              (rest vcf-files)))]
    (run-cmd out-file
             "zcat ~{(first vcf-files)} | ~{intersect-str} | "
             "vcfcreatemulti > ~{out-file}")
    (bgzip-index-vcf out-file)))
