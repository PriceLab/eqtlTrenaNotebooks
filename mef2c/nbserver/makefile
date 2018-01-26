build:
	docker build -t eqtl_mef2c_nbserver .

bash:
	docker run -ti --rm eqtl_mef2c_nbserver bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/jovyan/sharedData eqtl_mef2c_nbserver bash

run:
	docker run -d -p 8888:8888 eqtl_mef2c_nbserver start-notebook.sh --NotebookApp.token=""
	open http://localhost:8888

stop:
	- docker ps | grep eqtl_mef2c_nbserver | awk '{print $$1}' | xargs docker stop



