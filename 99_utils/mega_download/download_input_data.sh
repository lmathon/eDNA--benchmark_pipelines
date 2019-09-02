###############################################################################
## download data from MEGA server

## curl required
## wget required
## openssl required

echo "#########################################################################"
###############################################################################
## download reference data_base
URL="https://mega.nz/#"'!'"seonxSqJ"'!'"Alr84wbCCwJzRvMU25eMtYCKq11UOHTbHtFmm3k63xk"
FOLDER=00_Input_data/reference_database

if [[ -d "$FOLDER" ]]; then
    echo "$FOLDER already exists"
    echo "Remove this folder to download reference database"
else 
    echo "Downloading reference database from MEGA to "$FOLDER"..."
    bash 99_utils/mega_download/mega_fetch.sh $URL > 99_utils/mega_download/mega_info.txt
	FILE_URL=`sed '1q;d' 99_utils/mega_download/mega_info.txt`
	FILE_NAME=`sed '2q;d' 99_utils/mega_download/mega_info.txt`
	HEX=`sed '3q;d' 99_utils/mega_download/mega_info.txt`
	RAW_HEX=`sed '4q;d' 99_utils/mega_download/mega_info.txt`

	wget -O "${FILE_NAME}" "${FILE_URL}"
	cat "${FILE_NAME}" | openssl enc -d -aes-128-ctr -K "${HEX}" -iv "${RAW_HEX}" > "${FILE_NAME}".new
	mv -f "${FILE_NAME}".new "${FILE_NAME}"
	unzip "${FILE_NAME}"
	mv ${FILE_NAME/.zip/} 00_Input_data/
fi

echo "#########################################################################"
###############################################################################
## download fastq.gz files
URL="https://mega.nz/#"'!'"JXJSCAaa"'!'"uDWGamyHfEJkYt6BorxDetpu18EcRnmIyjzREJlOuXM"
FOLDER=00_Input_data/forward_reverse_reads

if [[ -d "$FOLDER" ]]; then
    echo "$FOLDER already exists"
    echo "Remove this folder to download grinder simulation reads"
else 
    echo "Downloading grinder simulations reads from MEGA to "$FOLDER"..."
	bash 99_utils/mega_download/mega_fetch.sh $URL > 99_utils/mega_download/mega_info.txt

	FILE_URL=`sed '1q;d' 99_utils/mega_download/mega_info.txt`
	FILE_NAME=`sed '2q;d' 99_utils/mega_download/mega_info.txt`
	HEX=`sed '3q;d' 99_utils/mega_download/mega_info.txt`
	RAW_HEX=`sed '4q;d' 99_utils/mega_download/mega_info.txt`

	wget -O "${FILE_NAME}" "${FILE_URL}"
	cat "${FILE_NAME}" | openssl enc -d -aes-128-ctr -K "${HEX}" -iv "${RAW_HEX}" > "${FILE_NAME}".new
	mv -f "${FILE_NAME}".new "${FILE_NAME}"
	unzip "${FILE_NAME}"
	mv ${FILE_NAME/.zip/} 00_Input_data/
fi