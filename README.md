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

### Detailed steps for UNIX

1. Create a new mySQL database called 'phyloisland'
2. Create a user 'pi' that has access to this database

```
sudo mysql -p
CREATE DATABASE phyloisland
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


## Usage

1. Upload FASTA files through the Upload tab
2. Annotate the genomes through Sequence Records tab
3. Any profiles built using the "Build profile" option are written to the Profiles tab

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D


## Credits

Gabe Foley and Joseph Box
