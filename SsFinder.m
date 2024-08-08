classdef SsFinder
    methods(Static)
        function best_match = findPss(samples,from_freq,to_freq,samples_per_symb)
            % finds primary SS in time-domain complex signal
            % returns structure of the match, containig next fields:
            % *NId2 - the part of NcellID (see 7.4.2.1 of TS38.211)
            % *corr - correlation between generated signal and `samples`
            % *lags - correlation `X` axis
            % *max_abs - maximum absolute value of the correlation func
            % *freq - subcarier index where starts PSS
            % *kSSB - subcarier shift to the SSB

            best_match=struct("max_abs",-1);
            
            % foreach NID2
            for id=0:2
                % generate PSS by NID2
                pss=[PssGenerator.generatePssByCellInfo(id) zeros(1,samples_per_symb-127)];
                % foreach freq. shift
                for freq=(from_freq:to_freq)+56
                    corr=xcorr(samples,ifft(circshift(pss,freq)));
                    % check if current NID2 and shift are better
                    if best_match.max_abs<max(abs(corr))
                        best_match.NId2=id;
                        best_match.freq=freq;
                        best_match.corr=corr/sqrt(sum(abs(samples).^2)*sum(abs(pss).^2));
                        best_match.max_abs=max(abs(corr));
                    end
                end
            end
            [~,best_match.lags]=xcorr(samples,pss);
            best_match.kSSB=best_match.freq-56;
        end
        function [NId1,max_corr] = checkSss(samples, pss_index, kssb, NId2, samples_per_symb)
            % finds NID1 by complex signal in time domain, 
            % time offset to PSS, and other parameters
            % choosing area, where to check SSS
            area=samples(pss_index+2*samples_per_symb:pss_index+3*samples_per_symb);
            max_corr=-1;
            
            % foreach NID1 calculate NcellID
            for id=(0:335)*3+NId2
                % generate signal by NcellID
                sss=[SssGenerator.generateSssByCellInfo(id) zeros(1,samples_per_symb-127)];
                corr=xcorr(area,ifft(circshift(sss,56+kssb)));
                current=max(abs(corr));
                % check if current ID is better
                if current>max_corr
                    max_corr=current;
                    NId1=floor(id/3);
                end
            end
        end

        function peaks_i=findPeaks(samples,level)
            % finds peaks by level
            % finds indexes of max values in continuous parts that exceed than a given level
            arguments
                samples
                level = 0.9
            end
            trs=(level<=samples);
            peaks_i=[];
            max_val=-1;
            max_i=-1;
            is_seq=false;
            for i=1:length(samples)
                if trs(i)
                    if ~is_seq
                        is_seq=true;
                        max_val=samples(i);
                        max_i=i;
                    else
                        if (samples(i))>max_val
                            max_i=i;
                            max_val=samples(i);
                        end
                    end
                else
                    if is_seq
                        is_seq=false;
                        peaks_i(end+1)=max_i;
                    end
                end
            end
        end
        function [NCellId,kSSB,t_index,cropped]=processSignalByPeakNo(...
                samples,...
                from_freq,...
                to_freq,...
                samples_per_symb,...
                peak_No,....
                peak_level)
            % finds PSS, checks SSS, extracts offsets and 
            % NcellID from time-domain complex signal samples
            match=SsFinder.findPss(samples,from_freq,to_freq,samples_per_symb);
            peaks=match.lags(SsFinder.findPeaks(abs(match.corr),peak_level));
            peaks=peaks(peaks>0); % to exclude cropped part pss
            NId1=SsFinder.checkSss(samples,peaks(peak_No),match.kSSB,match.NId2,samples_per_symb);
            NCellId=NId1*3+match.NId2;
            kSSB=match.kSSB;
            t_index=peaks(peak_No)+1;
            cropped=samples(t_index:end);
        end
    end
end