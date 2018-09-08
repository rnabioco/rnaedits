#! /usr/bin/env python

import sys
import os
import gzip
import argparse
from subprocess import call

def download_files(mdata, dl_prog, dl_prog_path, ssh_key):
  
  ascp_cmd = [
        dl_prog_path,
        "-QT",
        "-l",
        "200m",
        "-P",
        "33001",
        "-i",
        ssh_key]
  
  wget_cmd = [
    dl_prog_path
  ]
  
  if not os.path.exists("brainrest"):
    os.makedirs(os.path.join("brainrest", "untrimmed"))
  if not os.path.exists("hypothalamus"):
    os.makedirs("hypothalamus")
  if not os.path.exists("medulla"):
    os.makedirs("medulla")
    
  downloaded_fqs = []
  with open(mdata) as f:
      for line in f:
          if line.split()[0].startswith("sample_name"):
              continue
          fields = line.split("\t")
          
          fq1 = fields[8]
          fq2 = fields[9]
          fqs_url = [x.strip() for x in fields[19].split(";")]
          
          region = fields[1]
          if region == "Forebrain":
            region = os.path.join("brainrest", "untrimmed")
          else:
            region = region.lower()
          
          fq1 = os.path.join(region, fq1)
          fq2 = os.path.join(region, fq2)
          
          if dl_prog == 'ascp':
            fq1_cmd = ascp_cmd + ["era-fasp@" + fqs_url[0], fq1]
            fq2_cmd = ascp_cmd + ["era-fasp@" + fqs_url[1], fq2]
          elif dl_prog == 'wget':
            fqs_url = [x.replace("fasp", "ftp://ftp") for x in fqs_url]
            fq1_cmd = wget_cmd + ["-O", fq1, fqs_url[0]]
            fq2_cmd = wget_cmd + ["-O", fq2, fqs_url[1]]
          else:
            sys.exit("unknown download program {}".format(dl_prog))
            
          print("downloading fastq {}".format(fq1))
          call(fq1_cmd)
          
          print("downloading fastq {}".format(fq2))
          call(fq2_cmd)

def main():
    
    parser = argparse.ArgumentParser(description="""
    This script will download all of the raw data and 
    place the files into organized directories
    """)

    parser.add_argument('-i',
                        '--metadata',
                         help ='metadata tsv file (defaults to "docs/fastq_metadata.tsv")',
                       required = False,
                       default = os.path.join("..", "..", "docs",
                           "fastq_metadata.tsv"))
    parser.add_argument('-d',
                          '--downloader',
                          help ="""either 'ascp' or 'wget' """,
                       required = True)
    parser.add_argument('-p',
                          '--prog_path',
                          help ="""path to downloader binary, defaults to name of program, 
                          i.e. 'ascp', 'wget' '""",
                       required = False)
    parser.add_argument('-s',
                          '--ssh',
                          help ="""path to aspera openssh key file (i.e something like
                          ".aspera/connect/etc/asperaweb_id_dsa.openssh") """,
                       required = False, 
                       )
                      
    args=parser.parse_args()
    
    mdata = args.metadata
    dl_prog = args.downloader
    dl_prog_path = args.prog_path
    ssh_key = args.ssh
    
    if dl_prog == "ascp" and ssh_key is None:
      sys.exit("-s required for using aspera")
    if dl_prog_path is None:
      dl_prog_path = dl_prog
      
    download_files(mdata, dl_prog, dl_prog_path, ssh_key)

if __name__ == '__main__': main()
