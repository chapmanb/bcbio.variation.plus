(defproject bcbio.variation.plus "0.0.1-SNAPSHOT"
  :description "Extended version of bcbio.variation toolkit with additional utilities"
  :url "https://github.com/chapmanb/bcbio.variation.plus"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [bcbio.variation "0.1.3-SNAPSHOT"]
                 ;; GATK dependencies
                 [it.unimi.dsi/fastutil "6.5.3"]]
  :jvm-opts ["-Xms750m" "-Xmx2g"]
  :plugins [[lein-midje "3.1.3"]]
  :profiles {:dev {:dependencies
                   [[midje "1.6.0"]]}}
  :aot [bcbio.variation.plus]
  :main bcbio.variation.plus)
