## ImmunoPhenoDB

ImmunoPhenoDB is a Python library for creating and maintaining proteotrancriptomic (e.g. CITE-seq) data using a MySQL database. This library implements functionality for initializing the database schema, uploading new datasets, and general maintenance. This repository also contains all necessary backend infrastructure to access this data through APIs. Hosting an ImmunoPhenoDB server currently requires the usage of an instance through a cloud service provider such as AWS, GCP, Azure, etc.

## Cloud Instance Requirements
The instance must be running Ubuntu 22.04+ with Docker 26.0+ and MySQL 8.0+ installed. From the administrative console, ensure the following ports are opened:
1. Port 80 (HTTP)
2. Port 22 (SSH)
3. Port 8888 (JupyterLab)

## ImmunoPhenoDB Setup 
Loading datasets into the MySQL database is facilitated through the usage of a JupyterLab notebook. This notebook can be installed through the provided docker image, using the default port 8888:

```commandline
docker run -it -p 8888:8888 -v C:\Users\myusername\Documents\myfolder:/home/jovyan/work camaralab/python3:immunophenodb
```

Once the Docker container is established, connect to the JupyterLab notebook using a browser at ```http://<instance IP address>:8888```.


Create a folder called ```mysql_files``` in the home directory. Inside this directory, create a file called ```config.ini```, which will contain the credentials to connect to the MySQL database. An example ```config.ini``` file can be found in the repository.

Once the JupyterLab notebook is setup and configured to connect to the MySQL database, datasets must be normalized using instructions found in [Tutorial 2](https://immunopheno.readthedocs.io/en/main/notebooks/Example_2.html). 

Uploading a normalized dataset requires calling the ```load_csv_database()``` function, which automatically creates and populates all relevant MySQL tables.

## Backend API Setup and Connection

Accessing the populated MySQL tables requires all backend services (Flask, NodeJS, NGINX, Redis) to be established. These services enable users to query the database using the [ImmunoPheno](https://github.com/CamaraLab/ImmunoPheno) package.

After the MySQL database has been populated with normalized datasets, create a ```docker-compose.yml``` file in the home directory of the instance. This ```docker-compose.yml``` file orchestrates the services as docker containers and can be found in the repository. 

In the home directory containing the ```docker-compose.yml``` file, install and initialize these containers using the following command:
```commandline
docker compose up -d
```

Once these containers have been successfully built, connect to the server as shown in [Tutorial 1](https://immunopheno.readthedocs.io/en/main/notebooks/Example_1.html) using ```http://<instance IP address>``` inside the ```ImmunoPhenoDB_Connect``` class. 

