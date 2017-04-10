#!/bin/bash

bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed ARF2_1000_1_neg.bed -fo ARF2_1000_1_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed ARF2_1000_2_neg.bed -fo ARF2_1000_2_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed ARF2_1000_3_neg.bed -fo ARF2_1000_3_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed ARF2_1000_4_neg.bed -fo ARF2_1000_4_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed ARF2_1000_pos.bed -fo ARF2_1000_pos.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed MP_1000_1_neg.bed -fo MP_1000_1_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed MP_1000_2_neg.bed -fo MP_1000_2_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed MP_1000_3_neg.bed -fo MP_1000_3_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed MP_1000_4_neg.bed -fo MP_1000_4_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed MP_1000_pos.bed -fo MP_1000_pos.fas

exit 0
