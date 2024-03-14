#! /usr/bin/python
import os
import sys
import re
if __name__=="__main__":
    PS=sys.argv[1]
    com_1='find . -name "hw'+str(sys.argv[1])+'_*"| egrep -E "(JPG|pdf|png|jpg)"'
    f=os.popen(com_1)
    A=f.readlines()
    file_names=[]
    length=len(A)
    for i in range(1,length+1):
        search_pattern='hw'+str(sys.argv[1])+'_'+str(i)+'.(JPG|pdf|png|jpg)'
        for j in A:
            if re.search(search_pattern,j):
                if j.find('pdf') > -1:
                    file_names.append(j)
                else:
                    j_pdf=re.sub('(JPG|png|jpg)',"pdf",j)
                    f=os.popen('du '+j)
                    com_2='convert '+j[:-1]+' '+j_pdf
                    file_names.append(j_pdf) 
                    os.system(com_2)
                    com_4='rm '+j
                    os.system(com_4)
    com_3='pdftk '
    for file_name in file_names:
        com_3=com_3+file_name[:-1]+' '
    pdf_name='group5_PS'+sys.argv[1]+'_Xiehang_Yu.pdf'
    com_3=com_3+'cat output out.pdf' 
    os.system(com_3)
    os.system('cpdf -scale-to-fit a4portrait out.pdf -o '+pdf_name)
    os.system('rm out.pdf')
    for file_name in file_names:
        com_5='rm '+file_name
        os.system(com_5)
