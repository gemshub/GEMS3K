## kva2json

*kva2json* the CLI utility kva2json  to convert GEMS3K text I/O files from the old/present key-value arrays KVA format into the new JSON format (to be used e.g. in GEMS-Reaktoro) and back from JSON to key-value format.



## Application

```sh
Usage: kva2json [ option(s) ] -i|--import-from LST_IMPORT -e|--export-to LST_EXPORT [ -l|--dbr-list DBR_LIST ]
Recalculate task and export to other mode
Options:
	-h,	--help		show this help message

	-j,	--json        	write IPM, DCH and DBR files in json mode (default) 
	-t,	--key-value   	write IPM, DCH and DBR files in txt mode 
	-b,	--binary      	write IPM, DCH and DBR files in binary mode 

	-d,	--brife       	do not write data items that contain only default values (default false) 
	-c,	--comments    	write files with comments for all data entries ( in text mode ) (default false) 

```


#### Some run examples

* convert from key-value to json

```sh
> kva2json -j -i Kaolinite-in/pHtitr-dat.lst -e Kaolinite-json/pHtitr-dat.lst

```

* convert from json to key-value

```sh
> kva2json -t -c -i solvus-in/series1-dat.lst -e solvus-kv/series1-dat.lst

```

> Note: into very old *-ipm.dat* files need change
>  ```
>  # ID key of the initial chemical system definition
>  "Kaolinite G  pHtitr      0    0       1       25      0   "
>  ```
> to 
>  ```
>  # ID key of the initial chemical system definition
>  <ID_key> "Kaolinite G  pHtitr      0    0       1       25      0   "
> ```

