#!/usr/bin/env python3

import sqlite3 as sql

def main():
    os.system('echo > ./rnaseq.db')
    conn = sql.connect('./rnaseq.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE rnaseq_rmd_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path VARCHAR(255) NOT NULL
            );''')
    c.execute('''CREATE TABLE deseq2_output_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path VARCHAR(255) NOT NULL,
            rmd_file_id INT NOT NULL,
            FOREIGN KEY (rmd_file_id) REFERENCES rnaseq_rmd_files(file_id)
            );''')
    c.execute('''CREATE TABLE sig_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path VARCHAR(755) NOT NULL
            );''')
    conn.commit()

if __name__ == '__main__':
    main()

