input=
output=


obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $input | obisort -k count -r | obitab -o > $output