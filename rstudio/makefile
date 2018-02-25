build_trenaStudio:
	docker build -t trenastudio_01  -f Dockerfile-08  .

run:
	sudo docker run -d -p 8787:8787  -v /notebooks/:/home/rstudio/notebooks trenastudio_08


stop:
	docker stop `docker ps | tail -1 | awk '{print $$1}'`

buildBase:
	docker build -t pshannon/trenastudio .
	docker images | grep trenastudio


bash:
	docker run -ti --rm pshannon/trenastudio bash


packages:
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/BrowserViz-1)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trena)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trenaViz)


export:
	docker save -o trenastudio_08-image trenastudio_08:latest
	# scp trenastudio_08-image pshannon@newtrena:tmp/
	#newtrena.tmp> docker load < trenastudio_08-image
