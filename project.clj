(defproject bcbio.variation.plus "0.0.1-SNAPSHOT"
  :description "Extended version of bcbio.variation toolkit with additional utilities"
  :url "https://github.com/chapmanb/bcbio.variation.plus"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [bcbio.variation "0.1.1"]
                 ;; GATK dependencies
                 [it.unimi.dsi/fastutil "6.5.3"]
                 [org.clojure/core.incubator "0.1.3"]
                 [me.raynes/conch "0.5.0"]]
  :jvm-opts ["-Xms750m" "-Xmx2g"]
  :plugins [[lein-midje "3.0.1"]]
  :profiles {:dev {:dependencies
                   [[midje "1.5.1" :exclusions [org.clojure/clojure ordered]]]}}
  :aot [bcbio.variation.plus]
  :main bcbio.variation.plus)
