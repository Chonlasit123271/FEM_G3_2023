function [Coord,Quads,ELEMCon] = autoread(inputfile)
%%This function is an autoread of an inputfile
    % no of nodes is mentioned in 5th row and first column
    N_n      = dlmread(inputfile,'',[2-1 1-1 2-1 1-1]); %%Number of node
    N_e      = dlmread(inputfile,'',[4+N_n 0 4+N_n 0]); %%Number of element
    %%____NODE
    node_id     = dlmread(inputfile,'',[2 0 1+N_n 0]); %%Node id
    Coord1 = dlmread(inputfile,'',[2 0 1+N_n 3]); %%Node id with xyz coord
    Coord = dlmread(inputfile,'',[2 1 1+N_n 3]); %%xyz coord
    %%____ELEMENT
    ELEMCon1 = dlmread(inputfile,'',[5+N_n 0 5+N_n+N_e-1 12]);
    [row, col] = size(ELEMCon1);
    ELEMCon_stored =  dlmread(inputfile,'',[5+N_n 4 5+N_n+N_e-1 12]);
    for i = 1:size(ELEMCon_stored,1)
        if ELEMCon_stored(i,1) ==1
            columnnum_1(i) = i;
        elseif ELEMCon_stored(i,1) == 2
            columnnum_2(i-columnnum_1(length(columnnum_1))) = i;
        elseif ELEMCon_stored(i,1) == 4
            columnnum_4(i-columnnum_2(length(columnnum_2))) = i;
        elseif ELEMCon_stored(i,1) == 8
            columnnum_8(i-columnnum_4(length(columnnum_4))) = i;
        end
    end
    Quads = zeros(length(columnnum_4),4);
    ELEMCon = zeros(length(columnnum_8),8);
    for row = 1:length(columnnum_4)
        for col = 1:4
            Quads(row,col) = ELEMCon_stored(columnnum_4(row),col+1);
        end
    end
    for row = 1:length(columnnum_8)
        for col = 1:8
            ELEMCon(row,col) = ELEMCon_stored(columnnum_8(row),col+1);
        end
    end
end