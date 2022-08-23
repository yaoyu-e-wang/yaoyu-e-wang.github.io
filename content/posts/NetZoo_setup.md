---
title: "NetZooNotes"
date: 2022-06-24T15:08:42-04:00
draft: false
---

# Notes on setting up NetZoo Packages

## PyPanda

The only PyPanda working version is on Tian's GitHub account as the dev branch.  To set it up on an AWS EC2 instance running ubuntu 18.04, it needs to be install with python3 and it also loads 'setuptools' module that need to be installed.

An AWS machine is set up for testing with name NetZooTestVM with instance ID: **i-0fa07d648371d7adf**

Here is the code used for set up the environment:

```[sh]
sudo apt install python3
sudo apt install python3-pip
sudo python3-pip setuptools
```

Clone from Tian Wang's github page (https://github.com/twangxxx/pypanda-1) and run setup.

```[sh]
git clone https://github.com/twangxxx/pypanda-1.git
cd pypanda-1
python3 setup.py install

```

Use the following commands to run toy data:

```[sh]
# For PANDA Network
python 3 run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o output_panda.txt


# For Lioness Network
python3 run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o output_panda.txt -q output_lioness.txt

```
