build:
	docker build -t eqtl_mef2c_trena .
	docker images | grep eqtl_mef2c_trena

bash:
	docker run -ti --rm eqtl_mef2c_trena bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/trena/sharedData eqtl_mef2c_trena bash

packages:
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/BrowserViz-1)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trenaViz)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trena)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/eqtlTrenaNotebooks/mef2c/SingleGeneData)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/eqtlTrenaNotebooks/mef2c/MEF2C.data)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/eqtlTrenaNotebooks/mef2c/SingleGeneAnalyzer)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/MotifDb)

