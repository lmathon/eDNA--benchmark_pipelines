###############################################################################
## download data from MEGA server

## curl required
## wget required
## openssl required


## download reference data_base
URL="https://mega.nz/#"'!'"seonxSqJ"'!'"Alr84wbCCwJzRvMU25eMtYCKq11UOHTbHtFmm3k63xk"

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


## download fastq.gz files
URL="https://mega.nz/#"'!'"VaRGmAiK"'!'"Qx3e23heHBY3PcjclTGTFqg3KrHXVB2HdY_ZMhvQxak"

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