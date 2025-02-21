# How to run

The project itself is built only on C++ so the only thing needed to be done to compile and run the project (other than having the correct things installed) is run to run the following on the terminal, on the same directory as the file you are currently reading is in:

```bash
make
./bin/cmb
```

If you want to run the Python script to generate the plots, you can simply navegate to it and run it with:

```bash
cd src/python/
python plotting.py
```

If you want to run the Python script for specific calculations, you first need to compile the C++ code in order for Python to be able to access it too, by running the following on the terminal on the same directory as the file you are currently reading is in:

```bash
make python-lib
```

Then you can run it normally by again navigating to the script and then calling Python with:

```bash
cd src/python/
python calculations.py
```

Pro tip: If all else fails, try getting rid of everything first by running `make clean` and then try again.