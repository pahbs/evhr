#!/usr/bin/env python
#
# Query nga inventory on ADAPT by catalog ID
#
import psycopg2
import csv
import argparse
import os, errno
import shutil
import os

def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)

def query_db_catid_test(catID, prod_code='M1BS', out_dir='/explore/nobackup/people/pmontesa', db_table='nga_footprint_master_V2'):
    '''Query and select scenes from latest database
    '''
    with psycopg2.connect(database="arcgis", user="pmontesa", password=os.environ['NGADBPASS'], host="arcdb04", port="5432") as dbConnect:

        cur = dbConnect.cursor() # setup the cursor
        selquery =  "SELECT S_FILEPATH, SENSOR, CATALOG_ID, ACQ_TIME FROM %s WHERE CATALOG_ID = '%s' AND PROD_CODE = '%s'" %(db_table, catID, prod_code)
        #selquery =  "SELECT * FROM %s WHERE CATALOG_ID = '%s' AND PROD_CODE = '%s'" %(db_table, catID, prod_code)
        cur.execute(selquery)
        selected=cur.fetchall()

    return selected

def query_db_catid(catID: str, prod_code: str, out_dir: str, out_csv_fn: str, username: str, password: str, db_table='nga_footprint_master_V2', symlink=True):
    """
    Returns to stdout the ADAPT dir that has the imagery associated with the input catID
    If out_dir specified, copies all images and xmls associated with catid into out_dir as symbolic links
    """
    
    if out_csv_fn is not None:
        f = open( out_csv_fn, 'w')
    
    # Search ADAPT's NGA database for catID
    imglist=[]

    #with psycopg2.connect(database="ngadb01", user="anon", host="ngadb01", port="5432") as dbConnect:
    with psycopg2.connect(database="arcgis", user=username, password=password, host="arcdb04", port="5432") as dbConnect:

        cur = dbConnect.cursor() # setup the cursor

        #selquery =  "SELECT s_filepath, sensor, acq_time, cent_lat, cent_long FROM %s WHERE catalog_id = '%s' AND prod_code = '%s'" %(db_table, catID, prod_code)
        selquery =  "SELECT S_FILEPATH, SENSOR, ACQ_TIME FROM %s WHERE CATALOG_ID = '%s' AND PROD_CODE = '%s'" %(db_table, catID, prod_code)
        
        print( "\n\t Now executing database query on catID '%s' ..."%catID)
        cur.execute(selquery)
        selected=cur.fetchall()
        print( "\n\t Found '%s' scenes for catID '%s' \n"%(len(selected),catID))

        if len(selected) == 0:
            print('Exiting.')
            os._exit(1)

        if symlink:

            if not os.path.exists(out_dir): 
                 os.makedirs(out_dir)

            # This will only get the data that match the first prod_id, preventing replicated data from being copied. This should prevent mosaics from failing
            print(f'\t Making symlinks in output dir: {out_dir}')
            prod_id = selected[0][0].split('-')[-1].split('_')[0]
            print( "\t Creating symlinks for data associated with prod_id '%s'" %(prod_id))

            for i in range(0,len(selected)):

                if prod_id in selected[i][0]:
                    imglist.extend(selected[i][0])

                    #Copy ntf and xml to out_dir as symlinks
                    filename = os.path.split(selected[i][0])[1]
                    print(selected[i][0])

                    if out_csv_fn is not None:
                        f.write(selected[i][0]+'\n')
                    force_symlink( selected[i][0], os.path.join(out_dir, filename) )

                    try:
                        # shutil.Error: ... are the same file
                        # Just copy over the xmls, instead of creating a symlink to them
                        shutil.copy2(os.path.splitext(selected[i][0])[0]+'.xml', out_dir)
                    except Exception as e:
                        force_symlink( os.path.splitext(selected[i][0])[0]+'.xml', os.path.join(out_dir, os.path.splitext(filename)[0]+'.xml') )
            #return(imglist)

            # Print to stdout the dir from first record selected
            print(os.path.split(selected[0][0])[0])

    if out_csv_fn is not None:
        f.close()

def getparser():
    parser = argparse.ArgumentParser(description="Query NGAdb with a catalog id")
    parser.add_argument('catID', default=None, type=str, help='Input catid')
    parser.add_argument('-NGA_DB_USER', default=os.environ['USER'], type=str, help='Username for db arcgis on on host arcdb04')
    parser.add_argument('-NGA_DB_PASS', default=os.environ['NGADBPASS'], type=str, help='Password for db arcgis on on host arcdb04')
    parser.add_argument('-prod_code', default='P1BS', type=str, help='Image production code: P1BS or M1BS')
    parser.add_argument('-out_dir', default='/explore/nobackup/people/pmontesa', type=str, help='Output pairname dir for symlinks')
    parser.add_argument('-db_table', default='nga_footprint_master_v2', type=str, help='Specify the db table name in the database')
    parser.add_argument('-out_csv_fn', default=None, help='Output CSV of paths')
    parser.add_argument('--no-symlink', dest='symlink', action='store_false', help='Turn off symlinking into output dir')
    parser.set_defaults(symlink=True)
    return parser

def main():

    parser = getparser()
    args = parser.parse_args()

    catID = args.catID
    prod_code = args.prod_code
    out_dir = args.out_dir
    db_table = args.db_table
    out_csv_fn = args.out_csv_fn

    if args.NGA_DB_PASS is None:
        print('Needs password. Exiting.')
        os._exit(1)
        
    query_db_catid(catID, prod_code, out_dir, out_csv_fn, args.NGA_DB_USER, args.NGA_DB_PASS, db_table, args.symlink)

if __name__ == "__main__":
    main()
