# PRFuMS

高通量自动化组装流程

1. ### 软件依赖

   cd-hit-est   4.6.4

   https://github.com/weizhongli/cdhit

   fr-hit   0.7

   http://weizhong-lab.ucsd.edu/frhit/download.php

   velvet   1.2.10

   https://github.com/dzerbino/velvet

   phrap/cross_match    1.090518

   http://www.phrap.org/consed/consed.html#howToGet
   
   seqtk    1.3 (r106)
   
   https://github.com/lh3/seqtk
   
   Python模块：PyYaml

2. ### 数据格式要求

   原始双端文件：测序样本的序列已清除barcode，序列文件不分割，barcode标注在序列注释最后，用:分割
   
   @SeqID [describiton]:barcode
   
   可自行对fastqReader.py中的_read_one_record方法修改Barcode的识别方式

3. ### 配置文件说明

   config配置文件示例：config.yml
   格式为YAML格式，‘：’号后留一个空格

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
   # 是否输出无barcode匹配的reads
   mismatch: True
   # 一个样本的最大reads数，太大的数值会严重影响组装速度
   maxsize: 4000000
   # 运行时的最大进程数
   thread: 8
   # 自动删除运行过程中产生的临时文件
   auto_del: True
   # 日志记录等级，info = 0，debug = 1
   logger_level: 1
   # 压缩生成的序列文件，能节约80%左右的硬盘空间
   gzip: True
   ```
   
   
   barcode样例：barcode.txt
   
   ```
   ACTTAGAC  Sample_1
   ```


4. ### 使用说明

   ```
   python3 wrapper.py config_file step[0|1-6|N1:N2]
   
   步骤说明：
   Step 1： 序列文件分割，创建测序样本文件夹，FASTQ格式转换为FASTA，输出ID.fasta
   Step 2： 引物序列清除，输出ID.clean.fasta
   Step 3： 载体序列清除，输出ID.noVector.fasta
   Step 4： Velvet组装，输出ID.ForPhrap1.fasata
   Step 5： Phrap组装，输出ID.phrap.fasta
   
   0是运行全部步骤，2是只运行step2，2:4是从step2运行到step4.
   
   python3 wrapper.py config_1.txt 0
   
   python3 wrapper.py config_2.txt 2
   
   python3 wrapper.py config_3.txt 2:4
   ```
   
