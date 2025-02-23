function [seedPos,Peak] = SearchPeak(rhoMat,Vfa)
%SEARCHPEAK 此处显示有关此函数的摘要


detWay = DensityPeaksIR();

m = size(rhoMat, 1);
% 计算密度和距离
[delta] = iterationElection(detWay, Vfa);

% 查询密度峰的向量位置
[ classInitial ] = singularFind( detWay, Vfa, delta );
singularIndex = find( classInitial ~=  0 );
% 由向量位置转为矩阵位置
classCenterRows = mod( singularIndex, m );
classCenterRows(classCenterRows == 0) = m;
classCenterCols = ceil( singularIndex / m );
%
seedPos = [ classCenterCols, classCenterRows ];
% seedDelta= delta(singularIndex);

Peak.delta = delta;
Peak.singularIndex = singularIndex;
Peak.classCenterCols = classCenterCols;
Peak.classCenterRows = classCenterRows;

end

