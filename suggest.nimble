# Package
version     = "2.0.0"
author      = "Charles Blake"
description = "mmap-persistent spell checking algorithms"
license     = "MIT/ISC"
bin         = @["suggest", "tspell"]

# Deps
requires "nim >= 1.2.0", "cligen >= 1.10.0", "nio >= 0.7.12"
skipDirs = @["data"]
