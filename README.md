# Phylo Island

Phylo Island is a front-facing database that allows for curation and annotation of genomic sequences. It was specifically designed to interrogate large datasets consisting of multi-domain proteins within their genomic context.

## Dependencies for local deployment

1. Python 3
2. MySQL server
3. Python modules as described in requirements.txt

## Installation for local deployment

### Broad steps

1. Create a new mySQL database called 'phyloisland'
2. Create a user 'pi' that has access to this database
3. Create the database schema as described in install/biosqldb-mysql.sql
4. Install the python modules as described in requirements.txt
5. Run the phyloinit.py file to create the database
6. Run routes.py to start the application

### Detailed steps for UNIX

1. Create a new mySQL database called 'phyloisland'
2. Create a user 'pi' that has access to this database

```
sudo mysql -p
CREATE DATABASE phyloisland;
CREATE USER 'pi'@'localhost' IDENTIFIED BY 'password';
GRANT ALL on phyloisland.* to 'pi'@'localhost';
```

3. Create the database schema as described in install/biosqldb-mysql.sql

```
cd install
mysql -p phyloisland < biosqldb-mysql.sql
```

4. Install the python modules as described in requirements.txt
```
pip install -r requirements.txt
```

5. Run the phyloinit.py file to create the database

```
python phyloinit.py
```

6. Run the routes.py file to start the Flask application
```
python routes.py
```

7. Open a browser and navigate to http://127.0.0.1:7777/phyloisland_experimental

## Usage

1. Upload FASTA files through the Upload tab
2. Annotate the genomes through Sequence Records tab
3. Any profiles built using the "Build profile" option are written to the Profiles tab

## Creating multiple databases

You might want to create multiple databases to try out different configurations without removing all of your data.

The following steps show how to change to a new database named "newdatabase"

1. Replace the following line in phyloinit.py
```
server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="phyloisland")
```
with
```
server = BioSeqDatabase.open_database(driver="MySQLdb", user="pi", passwd="", host="localhost", db="newdatabase")
```
2. Replace the following line in servers.py
```
bio_server_name = "phyloisland"
```
with
```
bio_server_name = "newdatabase"
```
3. Follow the install steps but change "phyloisland" to something else, for example - "newdatabase"

Additionally, so that you can make changes to these files without them being read by git, enter the following in the command line so that git will skip over them. Do this before making changes to the files.

```
git update-index --skip-worktree phyloinit.py
git update-index --skip-worktree servers.py
```


## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D


## Credits

Gabe Foley and Joseph Box
