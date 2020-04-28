# PRFuMS

高通量自动化组装流程

1. ### 软件依赖

   cd-hit-est

   https://github.com/weizhongli/cdhit

   fr-hit

   http://weizhong-lab.ucsd.edu/frhit/download.php

   velvet

   https://github.com/dzerbino/velvet

   phrap/cross_match

   http://www.phrap.org/consed/consed.html#howToGet

2. ### 数据格式要求

   原始双端文件：测序样本的序列已清除barcode，序列文件不分割，barcode标注在序列注释最后面

   @SeqID [describiton]barcode

   

3. ### 配置文件说明

   config配置文件示例：config.txt

   ```
   # 双端文件，fastq格式，支持gz压缩
   FW_file: /home/vet/public/opt/data/PARFuMS/datasets/1.clean.fastq.gz
   RC_file: /home/vet/public/opt/data/PARFuMS/datasets/2.clean.fastq.gz
   # barcode文件
   BC_file: /home/vet/public/opt/data/PARFuMS/datasets/barcode.txt
   # 文件输出目录
   work_dir: /home/vet/public/opt/data/PARFuMS/output
   # 引物载体序列文件，fasta格式
   primer_file: /home/vet/public/opt/data/PARFuMS/datasets/primer.fasta
   vector_file: /home/vet/public/opt/data/PARFuMS/datasets/vector.fasta
   ```
   
   
   barcode样例：barcode.txt
   
   ```
   ACTTAGAC  Sample_1
   ```


4. ### 使用说明

   ```
   python3 wrapper.py config_file step[1:6|0|N1:N2]
   
   步骤说明：
   Step 1： 序列文件分割，创建测序样本文件夹，FASTQ格式转换为FASTA，输出ID.fasta
   Step 2： 引物序列清除，输出ID.clean.fasta
   Step 3： 载体序列清除，输出ID.noVector.fasta
   Step 4： Velvet组装，输出ID.ForPhrap1.fasata
   Step 5： Phrap组装，输出ID.phrap.fasta
   
   python3 wrapper.py config_1.txt 0
   
   python3 wrapper.py config_2.txt 2
   
   python3 wrapper.py config_3.txt 2:4
   ```
   
