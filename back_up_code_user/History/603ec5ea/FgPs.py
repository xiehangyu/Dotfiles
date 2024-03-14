import sys
import re
import os

command = "gem5.opt se.py --cmd=/home/xiehangyu/yxh/study/USTC/mm/computerarchtecturehomework/HW1/HW2/cs251a-microbench-master/mm --cpu-type=O3CPU --mem-type=DDR3_1600_8x8 --l1d_repl=NMRURP --l1d_assoc=2 --l1d_size=8kB --caches"

replacements = {
    "NMRURP": ["NMRURP", "RandomRP","LRURP","LIPRP"],  # 要替换的l1d_repl参数值
    "assoc": ["2", "4", "8"],  # 要替换的l1d_assoc参数值
}

commands = []
for repl1 in replacements["NMRURP"]:
    for repl2 in replacements["assoc"]:
        replaced_cmd = re.sub(r"--l1d_repl=\w+", f"--l1d_repl={repl1}", command)
        replaced_cmd = re.sub(r"--l1d_assoc=\d+", f"--l1d_assoc={repl2}", replaced_cmd)
        print(replaced_cmd)
        os.system(replaced_cmd)
        command2="cp -r m5out /home/xiehangyu/yxh/study/USTC/mm/computerarchtecturehomework/HW1/HW3_2/"+repl1+"_"+repl2+".csv"

        os.system(command2)
    


