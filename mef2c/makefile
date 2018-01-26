build:
	docker-compose build

up:
	docker-compose --no-ansi up

clean:
	- docker-compose rm --force --stop
	- docker images | grep eqtl_mef2c_trena | awk '{print $$3}' | xargs docker rmi --force

