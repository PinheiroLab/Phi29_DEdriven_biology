using FASTX
using BioSequences
using StatsBase
using DataFrames
using CSV
using Distributions

# Unique ID sequences that mark start and end of each indel loop library
exo_seq1 = dna"TTTGAAACCACC"
exo_seq2 = dna"ATTGAAGATCAC" 
#Amplicon: 45bp, 15amino acids

tpr2_seq1 = dna"GGTAAAGTTCCT"
tpr2_seq2 = dna"GAAACCAAAGAT"
#Amplicon: 45bp, 15amino acids

thumb_seq1 = dna"AGCCGTAAAATG"
thumb_seq2 = dna"GATACCTTTACG"
#Amplicon: 45bp, 15amino acids


# Loading FASTA file
function load_fasta_file(file_name::String)
    seqs = []
    reader = open(FASTA.Reader, file_name)
    for record in reader
        seq = FASTX.FASTA.sequence(record)
        push!(seqs, seq)
    end
    close(reader)
    return seqs
end

# Finding library region based on IDs
function find_and_trim(data, ID1, ID2)
    seqs = []
    for i in 1:length(data)
        pos1 = findfirst(ExactSearchQuery(ID1), data[i])
        pos2 = findfirst(ExactSearchQuery(ID2), data[i])
        if isnothing(pos1) || isnothing(pos2)
            continue
        end    
        trimmed_seq = data[i][pos1.stop+1:pos2.start-1]
        push!(seqs,trimmed_seq)
    end
        
    return seqs
end

# Transating to protein, ignoring sequences with frameshifts
function dna_2_aa(DNA_seqs)
    aa_seqs = []
    for i in 1:length(DNA_seqs)
        rna_seq = convert(LongRNA{2}, DNA_seqs[i])
        if length(rna_seq) % 3 == 0
            push!(aa_seqs,translate(rna_seq, allow_ambiguous_codons=true, alternative_start=false))    
        end
    end
    return aa_seqs
end

#=
Quantifying unique InDels and their enrichment through selection
    The following code calculates different scores to analyse the effectiveness of selection and identify significantly enriched variants with potentially 
    enhanced XNA synthesis activity.

1. Frequency: counts are divided by the total number of protein sequences
    freq = Counts_R0/total_counts

2. Enrichment: Enrichment of each protein sequence through selection
    enrichment = log(frequency_R1/frequency_R0)

3. Enrichment_complex: this score takes into consideration the frequency of sequences before and after selection as well their 
   joint counts to fish out the highly represented/abundant sequences 
    EC = (frequency_R1 - frequency_R0)x((counts_R0 + counts_R1)/2)

4. Statistical E-test: This adapts Andrew Leifer's MATLAB code to implement Krishnamoorty's E-test [Source: Krishnamoorthy, K and Thomson, J. (2004). 
   A more powerful test for comparing two Poisson means. Journal of Statistical Planning and Inference, 119, 249-267].
   Computes the p-value of the unconditional test for testingone and two-sided hypotheses about the means of two Poisson distributions.
        INPUT:
        iside = side of the test; true for right-sided, false for two-sided 
        alpha = nominal level of the test 
        ki = count of the ith population, i = 1,2 
        ni = sample size from the ith population, i=1,2 
        d = the difference mean1 - mean2 under the H0
        OUTPUT: 
        p-value = p-value of the unconditional test

=#

function indel_classification(aa_seqs,wt_length)
    output=[]
    total = length(aa_seqs)
    unique_counts = countmap(aa_seqs)
    for (keys,vals) in unique_counts
        size = length(keys) - wt_length
        freq = vals/total
        push!(output, (keys,size,vals,freq))
    end
    output = DataFrame([[output[k][kk] for k in 1:length(output)] for kk in 1:length(output[1])], [:Seq, :InDel, :Count, :Freq])
    return sort!(output, [:Freq],rev=true)
end

function enrichment_complex(R0_freq,R1_freq,R0_counts,R1_counts)
        EC = (R1_freq .- R0_freq).*((R0_counts.+R1_counts)/2)
        return EC
end

function etest(data::DataFrame, iside::Bool=false)
    pvalues = []
    for i in 1:size(data)[1]
        k1 = data[!,:Count_R0][i]
        k2 = data[!,:Count_R1][i]
        n1 = sum(data[!,:Count_R0])
        n2 = sum(data[!,:Count_R1])
        pvalue = testPoissonSignificance(k1,k2,n1,n2)
        push!(pvalues,pvalue)
    end
  return pvalues
end

function testPoissonSignificance(k1::Int64, k2::Int64, n1::Int64, n2::Int64, d::Float64=0.0, iside::Bool=false)
  elhatk = (k1+k2)/(n1+n2)-d*n1/(n1+n2)
    var = (k1/ (n1^2) + k2/(n2^2))
    t_k1k2 = (k1/n1-k2/n2-d)/sqrt(var)
  pvalue=poistest(n1, n2, elhatk, t_k1k2, d, iside)

  return pvalue
end

function poistest(n1::Int64, n2::Int64, elhatk, t_k1k2, d::Float64=0.0, iside::Bool=false)

  # computing estimates of el1*n1 and el2*n2 under H_0
  pvalue=0 
  elhat1=n1*(elhatk+d)
  elhat2 = n2*elhatk

  # computing the modes 
    i1mode = floor(elhat1)
    i2mode = floor(elhat2)

  # initializing the probability at the i1mode
    pi1mode = pdf(Poisson(elhat1), i1mode)
    pi1 = pi1mode

  # initializing the probability at the i2mode
    pi2mode = pdf(Poisson(elhat2), i2mode)
    for i1 = i1mode:1000  
        if pi1 >= 1e-7       
          pvalue=sumi2(n1, n2, elhat2, t_k1k2, i1, pi1, i2mode, pi2mode, d, pvalue, iside)
          pi1 = elhat1*pi1/(i1+1)
        end
    end 

  # Label #1 
  i1 = i1mode-1
  pi1 = pi1mode
    pi1 = i1mode*pi1/elhat1

    for i1 = i1mode-1:-1: 0
      if pi1 >= 1e-7 
          pvalue=sumi2(n1, n2, elhat2, t_k1k2, i1, pi1, i2mode, pi2mode, d, pvalue, iside)
          pi1 = i1*pi1/elhat1
      end
    end
  
  return pvalue
end

function sumi2(n1, n2, elhat2, t_k1k2, i1, pi1, i2mode, pi2mode, d, pvalue, iside)
  pvalue = pvalue
  pi2 = pi2mode

    for i2 = i2mode:1000
        if pi2 >= 1.0e-07
            elhati1 = 1.0e0*i1/n1
            elhati2 = 1.0e0*i2/n2
      
            diffi = elhati1 - elhati2 - d
            var = 1.0e0*elhati1/n1 + 1.0e0*elhati2/n2
    
        if iside == true     
        if 1.0e0*i1/n1 - 1.0e0*i2/n2 <= d
                t_i1i2 = 0.0e0
        else
                t_i1i2 = diffi/sqrt(var)
        end
            
        if t_i1i2 >= t_k1k2 
          pvalue += pi1*pi2
        end
          
      elseif iside == false
              if abs(1.0e0*i1/n1 - 1.0e0*i2/n2) <= d 
                t_i1i2 = 0.0e0
              else
                t_i1i2 = diffi/sqrt(var)
        end
            
        if abs(t_i1i2) >= abs(t_k1k2) 
          pvalue += pi1*pi2
        end
      end
    end
          
    pi2 = elhat2*pi2/(i2+1.0e0)
    end

  i2 = i2mode-1
        pi2 = pi2mode
        pi2 = i2mode*pi2/elhat2

        for i2 = i2mode-1:-1: 0
          if pi2 >= 1.0e-07 
            elhati1 = 1.0e0*i1/n1
            elhati2 = 1.0e0*i2/n2
            diffi = elhati1 - elhati2 - d
            var = (1.0e0*elhati1/n1 + 1.0e0*elhati2/n2)

            if iside == true 
              if 1.0e0*i1/n1 - 1.0e0*i2/n2 <= d 
                t_i1i2 = 0.0e0
              else
                t_i1i2 = diffi/sqrt(var)
              end 

              if t_i1i2 >= t_k1k2 
              pvalue += pi1*pi2
            end

          elseif iside == false     
              if abs(1.0e0*i1/n1 - 1.0e0*i2/n2) <= d 
                t_i1i2 = 0.0e0
              else
                t_i1i2 = diffi/sqrt(var)
              end 

            if abs(t_i1i2) >= abs(t_k1k2) 
              pvalue += pi1*pi2
            end
            end 
        end
          pi2 = i2*pi2/elhat2

        end

  return pvalue
end

function indel_enrichment(R0,R1, wt_length)
    #counting sequences
    r0 = indel_classification(R0, wt_length)
    r1 = indel_classification(R1,wt_length)
    rename!(r0,:Count => :Count_R0, :Freq => :Freq_R0, :InDel => :InDel_R0)
    rename!(r1,:Count => :Count_R1, :Freq => :Freq_R1, :InDel => :InDel_R1)
    
    #merging dataframes of r0 and r1 sequences
    merged = innerjoin(r0, r1, on = :Seq)
    
    #calculating enrichment and enrichment complex of sequences
    enrichment = log.(merged[!,:Freq_R1] ./ merged[!,:Freq_R0])
    ec = enrichment_complex(merged[!,:Freq_R0],merged[!,:Freq_R1],merged[!,:Count_R0],merged[!,:Count_R1])
    merged[!,:Enrichment] = enrichment
    merged[!,:EC] = ec

    #calculating e-test scores 
    etest_score = etest(merged)
    merged[!,:ETest] = etest_score
    return  sort!(merged,[:EC,:ETest],rev=true) ## sorting output by EC, this can be modified to any other metric'
end


#EX0
exo_R0 =load_fasta_file("01_Datasets\\02_Preprocessed\\EXO-R0.fasta")
exo_R1 =load_fasta_file("01_Datasets\\02_Preprocessed\\EXO-R1.fasta")
exo_R0aa = dna_2_aa(find_and_trim(exo_R0, exo_seq1,exo_seq2))
exo_R1aa = dna_2_aa(find_and_trim(exo_R1, exo_seq1,exo_seq2))
exo_InDels_R0R1 = indel_enrichment(exo_R0aa,exo_R1aa,15)

CSV.write("03_Output\\exo_InDels_R0R1.csv", exo_InDels_R0R1)


#TPR2
tpr2_R0 =load_fasta_file("01_Datasets\\02_Preprocessed\\TPR2-R0.fasta")
tpr2_R1 =load_fasta_file("01_Datasets\\02_Preprocessed\\TPR2-R1.fasta")
tpr2_R0aa = dna_2_aa(find_and_trim(tpr2_R0, tpr2_seq1,tpr2_seq2))
tpr2_R1aa = dna_2_aa(find_and_trim(tpr2_R1, tpr2_seq1,tpr2_seq2))
tpr2_InDels_R0R1 = indel_enrichment(tpr2_R0aa,tpr2_R1aa,15)

CSV.write("03_Output\\tpr2_InDels_R0R1.csv", tpr2_InDels_R0R1)


#Thumb
thumb_R0 =load_fasta_file("01_Datasets\\02_Preprocessed\\thumb-R0.fasta")
thumb_R1 =load_fasta_file("01_Datasets\\02_Preprocessed\\thumb-R1.fasta")
thumb_R0aa = dna_2_aa(find_and_trim(thumb_R0, thumb_seq1,thumb_seq2))
thumb_R1aa = dna_2_aa(find_and_trim(thumb_R1, thumb_seq1,thumb_seq2))
thumb_InDels_R0R1 = indel_enrichment(thumb_R0aa,thumb_R1aa,15)

CSV.write("03_Output\\thumb_InDels_R0R1.csv", thumb_InDels_R0R1)