# MCsquare

Fast Monte Carlo dose calculation algorithm for the simulation of PBS proton therapy.

## Compile the code (Linux)

The following instructions were tested on Ubuntu 20.04.2 LTS 64-bit.

**1. Install the Intel OneAPI compiler suite:**
To install the Intel toolkit, use these commands in your terminal:

```
sudo apt autoremove 'intel-*kit'  'intel-oneapi*'
cd /tmp
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-basekit intel-hpckit
```

More details available on: 
https://software.intel.com/content/www/us/en/develop/articles/installing-intel-oneapi-toolkits-via-apt.html

**2. download the source code of MCsquare:**
If git is not installed on your system yet, use these commands in your terminal:

```
sudo apt update
sudo apt install git
```

Then, clone the git repository on your computer:

```
git clone https://gitlab.com/openmcsquare/MCsquare.git
```

**3. Compile the code:**
Initialize the Intel compiler toolkit:
```
source /opt/intel/oneapi/setvars.sh
```

Move to the MCsquare folder and compile:
```
cd MCsquare
make all
```

## Run MCsquare (Linux)
You can test your executable by printing the version using this command in your terminal:
```
./MCsquare_linux -v
```
Configure your simulation in the config.txt file.
Then, run the MCsquare launcher from your terminal:
```
./MCsquare
```
The launcher will automatically call the executable that correspond to your computer hardware.

You can try to run MCsquare with the sample input data:
```
./MCsquare Sample_input_data/config.txt 
```



## Instructions for Windows
Install the Windows version of the Intel OneAPI toolkit (with the HPC module).

Compile MCsquare by running Makefile.bat
(You may need to adapt the path to librarary directories in Makefile.bat)

Configure the simulation in config.txt
Then, run MCsquare with the launcher: MCsquare.bat
