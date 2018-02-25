build_trenaStudio:
	docker build -t trenastudio_frd3 -f Dockerfile-frd3  .

run:
	sudo docker run -d -p 8787:8787  -v /notebooks/:/home/rstudio/notebooks trenastudio_frd3

stop:
	docker stop `docker ps | tail -1 | awk '{print $$1}'`

bash:
	docker run -ti --rm pshannon/trenastudio bash


packages:
	(cd ./pkgs; R CMD build --no-build-vignettes  ~/github/eqtlTrenaNotebooks/mef2c/SingleGeneData)
	(cd ./pkgs; R CMD build --no-build-vignettes  ~/github/eqtlTrenaNotebooks/mef2c/SingleGeneAnalyzer)
	(cd ./pkgs; R CMD build --no-build-vignettes  ~/github/eqtlTrenaNotebooks/mef2c/FRD3.data)
	(cd ./pkgs; R CMD build --no-build-vignettes  ~/github/trenaViz)


export:
	docker save -o trenastudio_frd3_image trenastudio_frd3:latest
	# scp trenastudio_08-image pshannon@newtrena:tmp/
	#newtrena.tmp> docker load < trenastudio_08-image
