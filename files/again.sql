CREATE TABLE profile_blob(
  uid SERIAL,
  name VARCHAR(255),
  mimetype VARCHAR(255),
  filename VARCHAR(255),
  profile BLOB,
  size INT(10),
  PRIMARY KEY (uid),
  UNIQUE (name));