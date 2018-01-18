# -*- coding: utf-8 -*-
"""
Created on Fri May  5 17:10:49 2017

@author: em1c14
"""
import numpy as np
import re,sys


header = 'java -jar /local/software/GATK/3.7/source/GenomeAnalysisTK.jar -T SelectVariants -R /temp/hgig/EXOME_DATA/REF_GENOMES/HG38/hs38.fa -V FINAL_GTYPED_BIALLELIC.vcf.gz' 
#header = 'java -jar /local/software/GATK/3.7/source/GenomeAnalysisTK.jar -T SelectVariants -R /temp/hgig/EXOME_DATA/REF_GENOMES/HG38/hs38.fa -V GT_NOD2.vcf.gz'
def list_provider(gene):
    ## Initialize the scripts for vcf selection
    #IBD
    out_ibd = open(gene+'_IBD.sh','w')
    out_ibd.write('#!/bin/sh\n')
    out_ibd.write('module load GATK/3.7\n')
    out_ibd.write(header + ' -L '+ gene + '.intervals')
    #Levin
    out_lev = open(gene+'_Lev.sh','w')
    out_lev.write('#!/bin/sh\n')
    out_lev.write('module load GATK/3.7\n')
    out_lev.write(header + ' -L '+ gene + '.intervals')

    #Full
    out_full = open(gene+'_FULL.sh','w')
    out_full.write('#!/bin/sh\n')
    out_full.write('module load GATK/3.7\n')
    out_full.write(header + ' -L '+ gene + '.intervals')
    

    ## Merge individuals per set using -sn 
    ibd = ' -sn PR0001 -sn PR0002 -sn PR0003 -sn PR0004 -sn PR0005 -sn PR0007 -sn PR0008 -sn PR0009 -sn PR0010 -sn PR0011 -sn PR0012 -sn PR0014 -sn PR0015 -sn PR0018 -sn PR0020 -sn PR0021 -sn PR0022 -sn PR0023 -sn PR0025 -sn PR0026 -sn PR0027 -sn PR0028 -sn PR0029 -sn PR0030 -sn PR0031 -sn PR0032 -sn PR0034 -sn PR0035 -sn PR0036 -sn PR0039 -sn PR0040 -sn PR0041 -sn PR0042 -sn PR0043 -sn PR0044 -sn PR0045 -sn PR0046 -sn PR0047 -sn PR0048 -sn PR0049 -sn PR0050 -sn PR0051 -sn PR0052 -sn PR0053 -sn PR0054 -sn PR0055 -sn PR0056 -sn PR0057 -sn PR0058 -sn PR0059 -sn PR0060 -sn PR0061 -sn PR0062 -sn PR0063 -sn PR0064 -sn PR0065 -sn PR0066 -sn PR0067 -sn PR0068 -sn PR0069 -sn PR0070 -sn PR0071 -sn PR0072 -sn PR0074 -sn PR0075 -sn PR0076 -sn PR0077 -sn PR0078 -sn PR0079 -sn PR0080 -sn PR0081 -sn PR0082 -sn PR0083 -sn PR0084 -sn PR0085 -sn PR0086 -sn PR0087 -sn PR0088 -sn PR0089 -sn PR0090 -sn PR0091 -sn PR0092 -sn PR0093 -sn PR0095 -sn PR0096 -sn PR0097 -sn PR0098 -sn PR0099 -sn PR0100 -sn PR0101 -sn PR0102 -sn PR0103 -sn PR0104 -sn PR0105 -sn PR0106 -sn PR0107 -sn PR0108 -sn PR0109 -sn PR0110 -sn PR0111 -sn PR0112 -sn PR0113 -sn PR0114 -sn PR0115 -sn PR0116 -sn PR0117 -sn PR0118 -sn PR0119 -sn PR0120 -sn PR0121 -sn PR0122 -sn PR0123 -sn PR0124 -sn PR0125 -sn PR0126 -sn PR0127 -sn PR0128 -sn PR0129 -sn PR0130 -sn PR0131 -sn PR0132 -sn PR0133 -sn PR0134 -sn PR0135 -sn PR0136 -sn PR0137 -sn PR0138 -sn PR0140 -sn PR0141 -sn PR0142 -sn PR0143 -sn PR0144 -sn PR0145 -sn PR0146 -sn PR0148 -sn PR0149 -sn PR0150 -sn PR0151 -sn PR0153 -sn PR0155 -sn PR0156 -sn PR0158 -sn PR0159 -sn PR0160 -sn PR0161 -sn PR0163 -sn PR0164 -sn PR0165 -sn PR0166 -sn PR0167 -sn PR0170 -sn PR0171 -sn PR0172 -sn PR0173 -sn PR0174 -sn PR0176 -sn PR0177 -sn PR0178 -sn PR0179 -sn PR0180 -sn PR0181 -sn PR0182 -sn PR0183 -sn PR0184 -sn PR0185 -sn PR0186 -sn PR0187 -sn PR0188 -sn PR0190 -sn PR0191 -sn PR0192 -sn PR0193 -sn PR0194 -sn PR0195 -sn PR0196 -sn PR0197 -sn PR0198 -sn PR0199 -sn PR0203 -sn PR0204 -sn PR0206 -sn PR0207 -sn PR0208 -sn PR0210 -sn PR0213 -sn PR0214 -sn PR0215 -sn PR0216 -sn PR0217 -sn PR0218 -sn PR0219 -sn PR0220 -sn PR0221 -sn SOPR0222 -sn SOPR0223 -sn SOPR0224 -sn SOPR0225 -sn SOPR0226 -sn SOPR0227 -sn SOPR0228 -sn SOPR0229 -sn SOPR0230 -sn SOPR0231 -sn SOPR0232 -sn SOPR0233 -sn SOPR0234 -sn SOPR0236 -sn SOPR0237 -sn SOPR0238 -sn SOPR0240 -sn SOPR0241 -sn SOPR0243 -sn SOPR0244 -sn SOPR0245 -sn SOPR0246 -sn SOPR0247 -sn SOPR0249 -sn SOPR0250 -sn SOPR0251 -sn SOPR0252 -sn SOPR0253 -sn SOPR0254 -sn SOPR0255 -sn SOPR0256 -sn SOPR0257 -sn SOPR0258 -sn SOPR0259 -sn SOPR0260 -sn SOPR0261 -sn SOPR0262 -sn SOPR0263 -sn SOPR0266 -sn SOPR0267 -sn SOPR0268 -sn SOPR0270 -sn SOPR0271 -sn SOPR0273 -sn SOPR0274 -sn SOPR0275 -sn SOPR0276 -sn SOPR0277 -sn SOPR0278 -sn SOPR0279 -sn SOPR0280 -sn SOPR0283 -sn SOPR0284 -sn SOPR0285 -sn SOPR0286 -sn SOPR0287 -sn SOPR0289 -sn SOPR0290 -sn SOPR0291 -sn SOPR0293 -sn SOPR0294 -sn SOPR0296 -sn SOPR0300 -sn SOPR0301 -sn SOPR0303 -sn SOPR0304 -sn SOPR0306 -sn SOPR0307 -sn SOPR0308 -sn SOPR0309 -sn SOPR0310 -sn SOPR0311 -sn SOPR0312 -sn SOPR0313 -sn SOPR0314 -sn SOPR0316 -sn SOPR0317 -sn SOPR0318 -sn SOPR0320 -sn SOPR0323 -sn SOPR0325 -sn SOPR0326 -sn SOPR0327 -sn SOPR0333 -sn SOPR0334 -sn SOPR0336 -sn SOPR0339 -sn SOPR0340 -sn SOPR0342 -sn SOPR0343 -sn SOPR0344 -sn SOPR0345 -sn SOPR0346 -sn SOPR0347 -sn SOPR0348 -sn SOPR0349 -sn SOPR0351 -sn SOPR0352 -sn SOPR0353 -sn SOPR0355 -sn SOPR0357 -sn SOPR0359 -sn SOPR0360 -sn SOPR0361 -sn SOPR0362 -sn SOPR0365 -sn SOPR0367 -sn SOPR0368 -sn SOPR0369 -sn SOPR0370 -sn SOPR0372 -sn SOPR0377 -sn SOPR0378 -sn SOPR0380 -sn SOPR0381 -sn SOPR0384'

    lev = ' -sn Ctrl_100 -sn Ctrl_101 -sn Ctrl_102 -sn Ctrl_103 -sn Ctrl_104 -sn Ctrl_105 -sn Ctrl_106 -sn Ctrl_107 -sn Ctrl_108 -sn Ctrl_109 -sn Ctrl_10 -sn Ctrl_110 -sn Ctrl_111 -sn Ctrl_112 -sn Ctrl_113 -sn Ctrl_114 -sn Ctrl_115 -sn Ctrl_116 -sn Ctrl_117 -sn Ctrl_118 -sn Ctrl_11 -sn Ctrl_120 -sn Ctrl_121 -sn Ctrl_122 -sn Ctrl_123 -sn Ctrl_124 -sn Ctrl_126 -sn Ctrl_127 -sn Ctrl_128 -sn Ctrl_129 -sn Ctrl_12 -sn Ctrl_130 -sn Ctrl_131 -sn Ctrl_132 -sn Ctrl_133 -sn Ctrl_134 -sn Ctrl_137 -sn Ctrl_138 -sn Ctrl_139 -sn Ctrl_13 -sn Ctrl_143 -sn Ctrl_144 -sn Ctrl_145 -sn Ctrl_148 -sn Ctrl_14 -sn Ctrl_150 -sn Ctrl_151 -sn Ctrl_152 -sn Ctrl_153 -sn Ctrl_154 -sn Ctrl_155 -sn Ctrl_156 -sn Ctrl_157 -sn Ctrl_158 -sn Ctrl_159 -sn Ctrl_15 -sn Ctrl_160 -sn Ctrl_162 -sn Ctrl_163 -sn Ctrl_164 -sn Ctrl_165 -sn Ctrl_166 -sn Ctrl_167 -sn Ctrl_168 -sn Ctrl_169 -sn Ctrl_16 -sn Ctrl_170 -sn Ctrl_171 -sn Ctrl_172 -sn Ctrl_174 -sn Ctrl_175 -sn Ctrl_177 -sn Ctrl_178 -sn Ctrl_179 -sn Ctrl_17 -sn Ctrl_181 -sn Ctrl_182 -sn Ctrl_183 -sn Ctrl_184 -sn Ctrl_185 -sn Ctrl_186 -sn Ctrl_187 -sn Ctrl_188 -sn Ctrl_189 -sn Ctrl_18 -sn Ctrl_190 -sn Ctrl_191 -sn Ctrl_192 -sn Ctrl_193 -sn Ctrl_194 -sn Ctrl_195 -sn Ctrl_196 -sn Ctrl_197 -sn Ctrl_198 -sn Ctrl_199 -sn Ctrl_19 -sn Ctrl_1 -sn Ctrl_200 -sn Ctrl_201 -sn Ctrl_202 -sn Ctrl_203 -sn Ctrl_204 -sn Ctrl_205 -sn Ctrl_206 -sn Ctrl_207 -sn Ctrl_209 -sn Ctrl_20 -sn Ctrl_211 -sn Ctrl_212 -sn Ctrl_213 -sn Ctrl_214 -sn Ctrl_215 -sn Ctrl_21 -sn Ctrl_22 -sn Ctrl_23 -sn Ctrl_24 -sn Ctrl_25 -sn Ctrl_26 -sn Ctrl_27 -sn Ctrl_28 -sn Ctrl_29 -sn Ctrl_2 -sn Ctrl_30 -sn Ctrl_31 -sn Ctrl_32 -sn Ctrl_33 -sn Ctrl_34 -sn Ctrl_35 -sn Ctrl_36 -sn Ctrl_37 -sn Ctrl_38 -sn Ctrl_39 -sn Ctrl_3 -sn Ctrl_40 -sn Ctrl_41 -sn Ctrl_42 -sn Ctrl_43 -sn Ctrl_44 -sn Ctrl_45 -sn Ctrl_46 -sn Ctrl_47 -sn Ctrl_48 -sn Ctrl_49 -sn Ctrl_4 -sn Ctrl_50 -sn Ctrl_51 -sn Ctrl_52 -sn Ctrl_53 -sn Ctrl_54 -sn Ctrl_55 -sn Ctrl_56 -sn Ctrl_57 -sn Ctrl_58 -sn Ctrl_59 -sn Ctrl_5 -sn Ctrl_60 -sn Ctrl_61 -sn Ctrl_62 -sn Ctrl_63 -sn Ctrl_64 -sn Ctrl_65 -sn Ctrl_66 -sn Ctrl_67 -sn Ctrl_68 -sn Ctrl_69 -sn Ctrl_6 -sn Ctrl_70 -sn Ctrl_71 -sn Ctrl_72 -sn Ctrl_73 -sn Ctrl_74 -sn Ctrl_75 -sn Ctrl_76 -sn Ctrl_77 -sn Ctrl_78 -sn Ctrl_79 -sn Ctrl_7 -sn Ctrl_80 -sn Ctrl_81 -sn Ctrl_82 -sn Ctrl_83 -sn Ctrl_84 -sn Ctrl_85 -sn Ctrl_86 -sn Ctrl_87 -sn Ctrl_88 -sn Ctrl_89 -sn Ctrl_8 -sn Ctrl_90 -sn Ctrl_91 -sn Ctrl_92 -sn Ctrl_93 -sn Ctrl_94 -sn Ctrl_95 -sn Ctrl_96 -sn Ctrl_97 -sn Ctrl_98 -sn Ctrl_99 -sn Ctrl_9'
    
    out_ibd.write(ibd)
    out_ibd.write(' -o ' + gene + '_IBD.vcf\n')
    
    out_lev.write(lev)
    out_lev.write(' -o ' + gene + '_Lev.vcf\n')
    
    
    out_full.write(ibd + lev)
    out_full.write(' -o ' + gene + '_FULL.vcf\n')
    
    out_ibd.close()
    out_lev.close()
    out_full.close()


list_provider(sys.argv[1])

