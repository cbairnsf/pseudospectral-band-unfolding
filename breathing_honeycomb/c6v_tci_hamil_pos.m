%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ssh_hamil_pos.m
% Alex Cerjan
% 2.10.22
%Descriptions of inputs added by CAB 08/29/25
%	nUCperSide	-- number of unit cells per side
%	t_2 		-- outside unit cell coupling
%	e_a, e_b	-- on-site energies
%	randSeed	-- seed for pseudorandom number generation
%	disStr		-- strength of disorder parameter
%Example usage: %	[HH, XX, YY, mirrorY, CC, mirrorX] = c6v_tci_hamil_pos(nUCperSide,t_2,e_a,e_b,randSeed,disStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function defn modified to remove irrelevant stuff CAB 08/30/25
%function [HH, XX, YY, mirrorY, CC, mirrorX] = c6v_tci_hamil_pos(nUCperSide,t_2,e_a,e_b,randSeed,disStr)
function [HH, XX, YY, T1, T2, T3] = c6v_tci_hamil_pos(nUCperSide, t_2)

	%Hardcoded values added by CAB 08/30/25
	e_a = 0;
	e_b = 0;

    % disorder parameters:
    %disStr = 0;
    %randSeed = 1;
    %rng(randSeed); % 'shuffle'
    
    %% lattice parameters:
    a = 3; % 2D lattice spacing
    % We're setting this to 3 so the site-site spacing within a unit cell = 1.
    
    %% initialization:
    coeffs = [t_2,e_a,e_b,a];

    %% total system setup:
    [HH] = genHamiltonian(nUCperSide,coeffs);

    [XX,YY] = genPosition(nUCperSide,coeffs);

	%Added by CAB b/c I need translation operators
	[T1, T2, T3] = genTranslation(XX, YY, a);

    %Everything below here irrelevant to CAB, so commented out 08/30/25
    %if (disStr > 0)
    %    assert(false);
    %    [disMat] = genDisorder(nUCperSide);
    %    HH = HH + disStr*disMat;
    %end

    %[CC] = genChiral(nUCperSide);
    %[mirrorY] = genMirrorY(nUCperSide);
    %[mirrorX] = genMirrorX(nUCperSide);
    %mirrorY = 0;

    %figure(1);
    %imagesc(HH);
    %figure(2);
    %imagesc(XX);
    %figure(3);
    %imagesc(Invers);
    %figure(4);
    %imagesc(Invers*XX*Invers + XX);

    %if (e_a == e_b)
    %    assert(isequal(HH*CC,-CC*HH));
    %    assert(isequal(HH*mirrorX,mirrorX*HH));
    %end
    %assert(isequal(XX*CC,CC*XX));
    %assert(isequal(YY*CC,CC*YY));
    %assert(isequal(XX*mirrorX,-mirrorX*XX));
    %assert(isequal(YY*mirrorX,mirrorX*YY));

    %assert(isequal(HH*mirrorY,mirrorY*HH));
    %assert(isequal(XX*mirrorY,mirrorY*XX));
    %assert(isequal(YY*mirrorY,-mirrorY*YY));

end

function [HH] = genHamiltonian(nUCperSide,coeffs)

    t_out = coeffs(1);
    e_a = coeffs(2);
    e_b = coeffs(3);
    %a = coeffs(4);

    t_in = 1;

	%commented out by CAB
   % if abs(t_out) > 1
   %     assert(t_out >= 0);
   %     t_in = 2-t_out;
   %     t_out = 1;
   % end
    
    % how many unit cells are there?
    nRows = 2*nUCperSide - 1;
    UCperRow = [nUCperSide:1:(2*nUCperSide-1),(2*nUCperSide-2):-1:nUCperSide];

    assert(length(UCperRow)==nRows);

    % how many non-zero values are there?
    nVals = 4*6*sum(UCperRow);
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvec = zeros(nVals,1);

    ii = 1;
    for rr=1:nRows
        for cc=1:UCperRow(rr)
            
            % we're going to start indexing the unit cells at 0,
            % that way, we can add to this for each lattice site.
            ucCur = sum(UCperRow(1:(rr-1))) + (cc-1);

            % The system is arrayed as:
            %        .   .   .   .
            %      .   .   .   .   .
            %    .   .   .   .   .   .
            %  .   .   .   .   .   .   .
            %    .   .   .   .   .   .
            %      .   .   .   .   .
            %        .   .   .   .  
            %  ---> x
            %  |
            %  v y

            % there are 6 potential other unit cells that might be adjacent to the current unit cell,
            % but we're only going to use 3, and put in both the forwards and backwards couplings simultaneously.
            ucPrevX = -1;
            if (cc>1)
                ucPrevX = ucCur-1;
            end
            %ucNextX = 0;
            %if (cc<UCperRow(rr))
            %    ucNextX = ucCur+1;
            %end

            ucPrevUL = -1;
            if (rr>1)
                if (rr <= ceil(nRows/2)) && (cc>1)
                    ucPrevUL = sum(UCperRow(1:(rr-2))) + (cc-2);
                end
                if (rr > ceil(nRows/2))
                    ucPrevUL = sum(UCperRow(1:(rr-2))) + (cc-1);
                end
            end
            ucPrevUR = -1;
            if (rr>1)
                if (rr <= ceil(nRows/2)) && (cc<UCperRow(rr))
                    ucPrevUR = sum(UCperRow(1:(rr-2))) + (cc-1);
                end
                if (rr > ceil(nRows/2))
                    ucPrevUR = sum(UCperRow(1:(rr-2))) + (cc-0);
                end
            end
            % there are 6 sites within each unit cell.
            % the numbering starts with site 1 at 3-o-clock, and continues clockwise. (i.e., site 2 is at 5-o-clock)
            idxCur1 = 6*ucCur + 1;
            idxCur2 = 6*ucCur + 2;
            idxCur3 = 6*ucCur + 3;
            idxCur4 = 6*ucCur + 4;
            idxCur5 = 6*ucCur + 5;
            idxCur6 = 6*ucCur + 6;

            idxPrevX1 = 0;
            if (ucPrevX ~= -1)
                idxPrevX1 = 6*ucPrevX + 1;
            end
            idxPrevUL2 = 0;
            if (ucPrevUL ~= -1)
                idxPrevUL2 = 6*ucPrevUL + 2;
            end
            idxPrevUR3 = 0;
            if (ucPrevUR ~= -1)
                idxPrevUR3 = 6*ucPrevUR + 3;
            end

            % on-site energies (A-sites are 1, 3, 5; B-sites are 2, 4, 6)
            xvec(ii) = idxCur1; yvec(ii) = idxCur1; vvec(ii) = e_a; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxCur2; vvec(ii) = e_b; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxCur3; vvec(ii) = e_a; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxCur4; vvec(ii) = e_b; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxCur5; vvec(ii) = e_a; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxCur6; vvec(ii) = e_b; ii = ii+1;

            % within unit cell couplings
            xvec(ii) = idxCur1; yvec(ii) = idxCur2; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxCur1; vvec(ii) = -t_in; ii = ii+1;

            xvec(ii) = idxCur2; yvec(ii) = idxCur3; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxCur2; vvec(ii) = -t_in; ii = ii+1;

            xvec(ii) = idxCur3; yvec(ii) = idxCur4; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxCur3; vvec(ii) = -t_in; ii = ii+1;

            xvec(ii) = idxCur4; yvec(ii) = idxCur5; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxCur4; vvec(ii) = -t_in; ii = ii+1;

            xvec(ii) = idxCur5; yvec(ii) = idxCur6; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxCur5; vvec(ii) = -t_in; ii = ii+1;

            xvec(ii) = idxCur6; yvec(ii) = idxCur1; vvec(ii) = -t_in; ii = ii+1;
            xvec(ii) = idxCur1; yvec(ii) = idxCur6; vvec(ii) = -t_in; ii = ii+1;

            % between unit cell couplings
            xvec(ii) = idxCur4; yvec(ii) = idxPrevX1; vvec(ii) = -t_out; ii = ii+1;
            xvec(ii) = idxPrevX1; yvec(ii) = idxCur4; vvec(ii) = -t_out; ii = ii+1;

            xvec(ii) = idxCur5; yvec(ii) = idxPrevUL2; vvec(ii) = -t_out; ii = ii+1;
            xvec(ii) = idxPrevUL2; yvec(ii) = idxCur5; vvec(ii) = -t_out; ii = ii+1;

            xvec(ii) = idxCur6; yvec(ii) = idxPrevUR3; vvec(ii) = -t_out; ii = ii+1;
            xvec(ii) = idxPrevUR3; yvec(ii) = idxCur6; vvec(ii) = -t_out; ii = ii+1;
        end % close cc loop
    end % close rr loop
    
    yvec(xvec == 0) = [];
    vvec(xvec == 0) = [];
    xvec(xvec == 0) = [];

    xvec(yvec == 0) = [];
    vvec(yvec == 0) = [];
    yvec(yvec == 0) = [];

    HH = sparse(xvec,yvec,vvec,6*sum(UCperRow),6*sum(UCperRow));

end % close function
    
function [disMat] = genDisorder(nUC)

    nVals = 2*2*nUC-2;
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvec = zeros(nVals,1);
    ii=1;

    for nn=1:nUC
        idxCurA = 2*nn-1;
        idxCurB = 2*nn;
        idxPrevB = 2*(nn-1);
        %idxNextA = 2*(nn+1)-1;

        disCur = (rand(1)-0.5);
        xvec(ii) = idxCurA; yvec(ii) = idxCurB; vvec(ii) = -disCur; ii = ii+1;
        xvec(ii) = idxCurB; yvec(ii) = idxCurA; vvec(ii) = -disCur; ii = ii+1;

        if nn > 1
            disCur = (rand(1)-0.5);
            xvec(ii) = idxCurA; yvec(ii) = idxPrevB; vvec(ii) = -disCur; ii = ii+1;
            xvec(ii) = idxPrevB; yvec(ii) = idxCurA; vvec(ii) = -disCur; ii = ii+1;
        end 
    end % close xx loop

    disMat = sparse(disX,disY,disV,2*nUC,2*nUC);
    
end

function [XX,YY] = genPosition(nUCperSide,coeffs)

    a = coeffs(4);

    nRows = 2*nUCperSide - 1;
    UCperRow = [nUCperSide:1:(2*nUCperSide-1),(2*nUCperSide-2):-1:nUCperSide];
    assert(length(UCperRow)==nRows);

    % how many non-zero values are there?
    nVals = 6*sum(UCperRow);
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvecX = zeros(nVals,1);
    vvecY = zeros(nVals,1);

    ii = 1;
    for rr=1:nRows
        for cc=1:UCperRow(rr)
            ucCur = sum(UCperRow(1:(rr-1))) + (cc-1);

            idxCur1 = 6*ucCur + 1;
            idxCur2 = 6*ucCur + 2;
            idxCur3 = 6*ucCur + 3;
            idxCur4 = 6*ucCur + 4;
            idxCur5 = 6*ucCur + 5;
            idxCur6 = 6*ucCur + 6;

            % there is always an odd number of unit cells in the "middle" row
            % and there are always an odd number of rows.
            % So, place origin in the center of the unit cell right at the lattice's center.

            ucX = a*((cc-1) - ((UCperRow(rr)-1)/2));
            ucY = (sqrt(3)/2)*a*((rr-1) - ((nRows-1)/2));

            siteSiteSpace = 1.05 * (a/3);
            
            xvec(ii) = idxCur1; yvec(ii) = idxCur1; vvecX(ii) = ucX + siteSiteSpace; vvecY(ii) = ucY; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxCur2; vvecX(ii) = ucX + (1/2)*siteSiteSpace; vvecY(ii) = ucY + (sqrt(3)/2)*siteSiteSpace; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxCur3; vvecX(ii) = ucX - (1/2)*siteSiteSpace; vvecY(ii) = ucY + (sqrt(3)/2)*siteSiteSpace; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxCur4; vvecX(ii) = ucX - siteSiteSpace; vvecY(ii) = ucY; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxCur5; vvecX(ii) = ucX - (1/2)*siteSiteSpace; vvecY(ii) = ucY - (sqrt(3)/2)*siteSiteSpace; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxCur6; vvecX(ii) = ucX + (1/2)*siteSiteSpace; vvecY(ii) = ucY - (sqrt(3)/2)*siteSiteSpace; ii = ii+1;

        end
    end

    vvecY = -1*vvecY;

    XX = sparse(xvec,yvec,vvecX,6*sum(UCperRow),6*sum(UCperRow));
    YY = sparse(xvec,yvec,vvecY,6*sum(UCperRow),6*sum(UCperRow));
end

function [CC] = genChiral(nUCperSide)

    nRows = 2*nUCperSide - 1;
    UCperRow = [nUCperSide:1:(2*nUCperSide-1),(2*nUCperSide-2):-1:nUCperSide];
    assert(length(UCperRow)==nRows);

    % how many non-zero values are there?
    nVals = 6*sum(UCperRow);
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvec = zeros(nVals,1);

    ii = 1;
    for rr=1:nRows
        for cc=1:UCperRow(rr)
            ucCur = sum(UCperRow(1:(rr-1))) + (cc-1);

            idxCur1 = 6*ucCur + 1;
            idxCur2 = 6*ucCur + 2;
            idxCur3 = 6*ucCur + 3;
            idxCur4 = 6*ucCur + 4;
            idxCur5 = 6*ucCur + 5;
            idxCur6 = 6*ucCur + 6;

            xvec(ii) = idxCur1; yvec(ii) = idxCur1; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxCur2; vvec(ii) = -1; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxCur3; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxCur4; vvec(ii) = -1; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxCur5; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxCur6; vvec(ii) = -1; ii = ii+1;
        end
    end

    CC = sparse(xvec,yvec,vvec,6*sum(UCperRow),6*sum(UCperRow));
end

function [mirrorY] = genMirrorY(nUCperSide)

    nRows = 2*nUCperSide - 1;
    UCperRow = [nUCperSide:1:(2*nUCperSide-1),(2*nUCperSide-2):-1:nUCperSide];
    assert(length(UCperRow)==nRows);

    % how many non-zero values are there?
    nVals = 6*sum(UCperRow);
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvec = zeros(nVals,1);

    ii = 1;
    for rr=1:nRows
        for cc=1:UCperRow(rr)
            ucCur = sum(UCperRow(1:(rr-1))) + (cc-1);
            rrOther = nRows - (rr-1);
            ucOther = sum(UCperRow(1:(rrOther-1))) + (cc-1);

            idxCur1 = 6*ucCur + 1;
            idxCur2 = 6*ucCur + 2;
            idxCur3 = 6*ucCur + 3;
            idxCur4 = 6*ucCur + 4;
            idxCur5 = 6*ucCur + 5;
            idxCur6 = 6*ucCur + 6;

            idxOther1 = 6*ucOther + 1;
            idxOther2 = 6*ucOther + 2;
            idxOther3 = 6*ucOther + 3;
            idxOther4 = 6*ucOther + 4;
            idxOther5 = 6*ucOther + 5;
            idxOther6 = 6*ucOther + 6;

            xvec(ii) = idxCur1; yvec(ii) = idxOther1; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther1; yvec(ii) = idxCur1; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxOther6; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther6; yvec(ii) = idxCur2; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxOther5; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther5; yvec(ii) = idxCur3; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxOther4; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther4; yvec(ii) = idxCur4; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxOther3; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther3; yvec(ii) = idxCur5; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxOther2; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther2; yvec(ii) = idxCur6; vvec(ii) = 1; ii = ii+1;
        end
    end

    mirrorY = sparse(xvec,yvec,vvec,6*sum(UCperRow),6*sum(UCperRow));
end

function [mirrorX] = genMirrorX(nUCperSide)

    nRows = 2*nUCperSide - 1;
    UCperRow = [nUCperSide:1:(2*nUCperSide-1),(2*nUCperSide-2):-1:nUCperSide];
    assert(length(UCperRow)==nRows);

    % how many non-zero values are there?
    nVals = 6*sum(UCperRow);
    xvec = zeros(nVals,1);
    yvec = zeros(nVals,1);
    vvec = zeros(nVals,1);

    ii = 1;
    for rr=1:nRows
        for cc=1:UCperRow(rr)
            ucCur = sum(UCperRow(1:(rr-1))) + (cc-1);
            ucOther = sum(UCperRow(1:(rr-1))) + (UCperRow(rr)-1 -(cc-1));

            idxCur1 = 6*ucCur + 1;
            idxCur2 = 6*ucCur + 2;
            idxCur3 = 6*ucCur + 3;
            idxCur4 = 6*ucCur + 4;
            idxCur5 = 6*ucCur + 5;
            idxCur6 = 6*ucCur + 6;

            idxOther1 = 6*ucOther + 1;
            idxOther2 = 6*ucOther + 2;
            idxOther3 = 6*ucOther + 3;
            idxOther4 = 6*ucOther + 4;
            idxOther5 = 6*ucOther + 5;
            idxOther6 = 6*ucOther + 6;

            xvec(ii) = idxCur1; yvec(ii) = idxOther4; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther1; yvec(ii) = idxCur1; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur2; yvec(ii) = idxOther3; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther6; yvec(ii) = idxCur2; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur3; yvec(ii) = idxOther2; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther5; yvec(ii) = idxCur3; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur4; yvec(ii) = idxOther1; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther4; yvec(ii) = idxCur4; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur5; yvec(ii) = idxOther6; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther3; yvec(ii) = idxCur5; vvec(ii) = 1; ii = ii+1;
            xvec(ii) = idxCur6; yvec(ii) = idxOther5; vvec(ii) = 1; ii = ii+1;
            %xvec(ii) = idxOther2; yvec(ii) = idxCur6; vvec(ii) = 1; ii = ii+1;
        end
    end

    mirrorX = sparse(xvec,yvec,vvec,6*sum(UCperRow),6*sum(UCperRow));
end

%Function added by CAB 08/30/25 to create translation operators from XX and YY
function [T1, T2, T3] = genTranslation(X, Y, a)
	a1 = -[3/2;  sqrt(3)/2];	%First vector balanced triangular approach
	a2 =  [3/2; -sqrt(3)/2];	%Second vector
	tol = 1e-2; 			%Numerical tolerance 
	T1_list_rows = [ ];		%Start with nothing
	T1_list_columns = [ ];
	T2_list_rows = [ ];		%Start with nothing
	T2_list_columns = [ ];
	for index = 1:size(X, 1)
		for jndex = 1:size(X, 2)
			%This is giving the wrong results. Way too many connections!
			displacement = [X(index, index) - X(jndex, jndex); Y(index, index) - Y(jndex, jndex)];
			a1_test = (displacement - a1)' * (displacement - a1); 
			a2_test = (displacement - a2)' * (displacement - a2); 
			if a1_test < tol
				T1_list_columns = [T1_list_columns, index];
				T1_list_rows = [T1_list_rows, jndex];
			elseif a2_test < tol
				T2_list_columns = [T2_list_columns, index];
				T2_list_rows = [T2_list_rows, jndex];
			end
		end
	end
	T1 = sparse(T1_list_columns, T1_list_rows, ones(1, length(T1_list_columns)), size(X, 1), size(X, 2));
	T2 = sparse(T2_list_columns, T2_list_rows, ones(1, length(T2_list_columns)), size(Y, 1), size(Y, 2));
	T3 = T2.' * T1.' + T1.' * T2.';	%Include both "left then right" and "right then left" paths so we don't lose connections because we can't jump off the flake 
	T3(T3 > 1) = 1;			%Get rid of all the 2's introduced by the sum above. We're just correcting for the left/right edges of the flake. 
end
