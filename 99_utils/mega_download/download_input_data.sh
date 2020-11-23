###############################################################################
## download data from MEGA server

## curl required
## wget required
## openssl required

echo "#########################################################################"
###############################################################################
## download reference data_base
URL="https://mega.nz/#"'!'"seonxSqJ"'!'"Alr84wbCCwJzRvMU25eMtYCKq11UOHTbHtFmm3k63xk"
FOLDER=benchmark_simulated_dataset/00_Input_data/reference_database



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
mv ${FILE_NAME/.zip/} benchmark_simulated_dataset/00_Input_data/reference_database/


echo "#########################################################################"
###############################################################################
## download fastq.gz files
#URL="https://mega.nz/#"'!'"JXJSCAaa"'!'"uDWGamyHfEJkYt6BorxDetpu18EcRnmIyjzREJlOuXM"
URL="https://mega.nz/#"'!'"AS4XASKL"'!'"dTkQOMgfrOPmEIjzMEQgzB7YO5_M7uD5E3lLSjdNmGw"

FOLDER=benchmark_simulated_dataset/00_Input_data/forward_reverse_reads


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
mv ${FILE_NAME/.zip/} benchmark_simulated_dataset/00_Input_data/forward_reverse_reads
