"""
version 1.0
author: Felix Plasser
description: Procedures for reading and writing files.
"""

import os

class dict_plus(dict):
    """
    Extension of a dictionary where data can be read from an input file.
    """
    def __init__(self, file_name='', revert=False):
        if not file_name == '':
            if not revert:
                self.read_from_file(file_name)
            else:
                self.revert_read_from_file(file_name)
    
    def read_from_file(self,file_name):
        """
        Read data from a file for the dictionary.
        """
        r_file = open(file_name, 'r')
        for line in r_file:
            spaceind = line.find(' ')
            self[line[:spaceind]]=line[spaceind+1:-1]
            
    def revert_read_from_file(self,file_name):
        """
        Read data from a file for the dictionary. Keys and values are switched
        """
        r_file = open(file_name, 'r')
        for line in r_file:
            spaceind = line.find(' ')
            self[line[spaceind+1:-1]]=line[:spaceind]

def chmkdirs(dir):
    """
    Tries to change to the directory *dir* and creates it if it does not exist.
    """
    try:
        os.chdir(dir)
    except:
        os.makedirs(dir)
        os.chdir(dir)
        
def line_to_words(line, sep=[' '], convert=None):
    """
    For parsing files.
    Returns a list of words that were given in <line>, separated by <sep>.
    If line ends in'\n', this is rejected.
    If <not convert==None> values are converted according to convert
        e.g. <convert=float> for conversion into floating point numbers.
    """
    if line[-1] == '\n':
        line = line[:-1]
    ret_list = []
    tmp_str = ''
    for let in line:
        if not let in sep:
            tmp_str += let
        else:
            if not tmp_str == '':
                if not convert==None:
                    tmp_str = eval(convert+'(tmp_str)')
                ret_list.append(tmp_str)
            tmp_str = ''
    if not tmp_str == '':
        if not convert==None:
            tmp_str = eval(convert+'(tmp_str)')
        ret_list.append(tmp_str)
            
    return ret_list

def change_file_lines(file_name, ind_cont=[]):
    """
    Change lines in a file (lines are overwritten).
    <ind_cont> is a list with indices and the contents for this line, e.g. ind_cont=[[3,'test']] writes 'test' into line 3.
    """
    lines = open(file_name, 'r').readlines()
    for ind, cont in ind_cont:
        lines[ind-1]=cont + '\n'

    w_file = open(file_name, 'w')
    w_file.writelines(lines)
    w_file.close()
    

class table_maker:
    """
    Class for writing output in columns.
    """
    def __init__(self, col_widths, cut=True, replace_list=[]):
        """
        Enter the widths of the columns in list <col_widths>.
        <replace_list> contains a double, list of items to be replaced,
            e.g. replace_list=[['.',',']] for European decimal notation.
        If <cut==True> content is cut to fit into the columns.
        """
        self.col_widths = col_widths
        self.ret_string = '' # string to be returned
        self.cut = cut
        self.replace_list = replace_list

    def write_header_line(self, words):
        """
        Writes a line with list <words>, prepending a # sign
        """
        plus_string = '#'
        for i, word in enumerate(words):
            if self.cut:
                plus_string += str(word)[:(self.col_widths[i]-1)].ljust(self.col_widths[i])
            else:
                plus_string += str(word).ljust(self.col_widths[i])

        for old,new in self.replace_list:
            plus_string = plus_string.replace(old, new)
    
        self.ret_string += plus_string + '\n'

    def write_line(self, words):
        """
        Writes a line with list <words>.
        """
        plus_string = ' '
        for i, word in enumerate(words):
            if self.cut:
                plus_string += str(word)[:(self.col_widths[i]-1)].ljust(self.col_widths[i])
            else:
                plus_string += str(word).ljust(self.col_widths[i])

        for old,new in self.replace_list:
            plus_string = plus_string.replace(old, new)
    
        self.ret_string += plus_string + '\n'

    def return_table(self):
        """
        Return the table that has been written with write_line.
        """
        return self.ret_string

    def write_to_file(self, file_name):
        """
        Write the table to a file.
        """
        write_to_file(self.return_table(), file_name)
        
class csv_maker:
    """
    Class for mading a csv file.
    """
    def __init__(self, sep=',', replace_list=[]):
        """
        <sep> is the separator in the output file.
        <replace_list> contains a double, list of items to be replaced,
            e.g. replace_list=[['.',',']] for European decimal notation.
        """
        self.ret_string = '' # string to be returned
        self.sep = sep
        self.replace_list = replace_list

        self.line_start = True    # a line is just starting

    def write_word(self, word):
        """
        Add a new entry to the csv.
        """
        if not self.line_start:
            self.ret_string += self.sep + word
        else:
            self.line_start = False
            self.ret_string += word

    def new_line(self):
        """
        Start a new line.
        """
        self.ret_string += '\n'
        self.line_start = True

    def write_line(self, words):
        """
        Writes a line with list <words>.
        """
        for word in words:
            self.write_word(word)
        self.new_line()

    def return_csv(self):
        """
        Return the table that has been written with write_line.
        """
        return self.ret_string

    def write_to_file(self, file_name):
        """
        Write the table to a file.
        """
        write_to_file(self.return_table(), file_name)
                      
def write_to_file(string, file_name):
    """
    Write <string> to file with <file_name>.
    """
    w_file = open(file_name, 'w')
    w_file.write(string)
    w_file.close()
        
if __name__ == '__main__':
    cm = csv_maker(sep=';')
    cm.write_line(['a', 'bea', 'da'])
    cm.write_word('3')
    cm.write_word('r')
    cm.new_line()
    cm.write_word('t')
    print cm.return_csv()

    #cm.write_to_file('test.txt')

