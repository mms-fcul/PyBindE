#!/usr/bin/python3

import re

def read_nb(file):
  NB=[]
  lines_no_comments=[]
  clean_lines=[]
  with open(file) as nb:
      next(nb)
      next(nb)
      lines = nb.read().splitlines()
      
      for line in lines:
          clean_line = line.split(";")
          #print(fields)
          lines_no_comments.append(clean_line[0])
      for line in lines_no_comments:
          cleanerline= re.sub('\s+',' ',line)
          clean_lines.append(cleanerline.strip())
      #print(clean_lines)
      clean_lines = list(filter(None, clean_lines))
      for line in clean_lines:
          fields = line.split(' ')
          #print(fields)
          nb_dict = {
              'i':    fields[0],
              'j':    fields[1],
              'c6':   float(fields[3]),
              'c12':  float(fields[4])
            }
          NB.append(nb_dict)
  return NB

read_nb('nb.itp')

def pair_search(criteria, string):
  list = []
  for pair in NB:
      if pair[criteria] == string:
          list.append(pair)
  return list
#pair_search('i', "CH3")