--Warning: Needs a lot of refactoring. FarthestPointFromSurface unnecessary as function only used once. Simplify if statements. Use indexing methods that avoid 2D arrays.
--Old code from 2018











function DistFromSurface(normalVec, posOnSurface, posn) --planar surface, negative output dist -> behind surface
	return normalVec.Unit:Dot(posn - posOnSurface)
end

function DistFromLine(posOnLine1, posOnLine2, posn)
	return ((posn - posOnLine1):Cross(posn - posOnLine2)).Magnitude/(posOnLine2 - posOnLine1).Magnitude
end

function ConvexHull3D(setOfPoints) --Works based on Quick-Hull Algorithm
	if #setOfPoints < 4 then error('3D Convex Hull requires at least 4 points.') end
	
	local function FarthestPointFromSurface(p1, p2, p3, pointSet) --p1, p2, p3 all on same plane (beware this function returning either of p1, p2, p3 as furthest!)
		local farthestDistance, farthestKey = 0
		local normalVec, posOnSurface = (p2 - p1):Cross(p2 - p3), p1
		
		for i = 1, #pointSet do
			local testDist = DistFromSurface(normalVec, posOnSurface, setOfPoints[pointSet[i]])
			if testDist > farthestDistance then
				farthestDistance, farthestKey = testDist, i
			end
		end
		return pointSet[farthestKey], setOfPoints[pointSet[farthestKey]], farthestDistance --returns farthest key in pointset, point and dist
	end
	
	----Initial phase
	
	--Create initial simplex (tetrahedron)
	local minXkey, minYkey, minZkey, maxXkey, maxYkey, maxZkey = 1, 1, 1, 1, 1, 1
	for key, point in pairs(setOfPoints) do
		if point.X < setOfPoints[minXkey].X then
			minXkey = key
		end
		if point.Y < setOfPoints[minYkey].Y then
			minYkey = key
		end
		if point.Z < setOfPoints[minZkey].Z then
			minZkey = key
		end
		if point.X > setOfPoints[maxXkey].X then
			maxXkey = key
		end
		if point.Y > setOfPoints[maxYkey].Y then
			maxYkey = key
		end
		if point.Z > setOfPoints[maxZkey].Z then
			maxZkey = key
		end
	end
	
	local extremePoints = {minXkey, minYkey, minZkey, maxXkey, maxYkey, maxZkey}
	local baseline_maxdist, baseline_k1, baseline_k2 = 0
	for i = 1, 6 do
		for j = 1, 6 do
			local testDist = (setOfPoints[extremePoints[i]] - setOfPoints[extremePoints[j]]).Magnitude
			if testDist > baseline_maxdist then
				baseline_maxdist, baseline_k1, baseline_k2 = testDist, extremePoints[i], extremePoints[j]
			end
		end
	end
	local baseline_p1, baseline_p2 = setOfPoints[baseline_k1], setOfPoints[baseline_k2]
	local apexPointDist, apexKey = 0
	for i = 1, #setOfPoints do
		local testDist = DistFromLine(baseline_p1, baseline_p2, setOfPoints[i])
		if testDist > apexPointDist then
			apexPointDist, apexKey = testDist, i
		end
	end
	local apexPoint = setOfPoints[apexKey]
	local I_farthestDistance, tipKey = 0
	local I_normalVec, I_posOnSurface = (baseline_p2 - baseline_p1):Cross(baseline_p2 - apexPoint), baseline_p1
	
	for i = 1, #setOfPoints do
		local testDist = math.abs(DistFromSurface(I_normalVec, I_posOnSurface, setOfPoints[i]))
		if testDist > I_farthestDistance then
			I_farthestDistance, tipKey = testDist, i
		end
	end
	local tipPoint = setOfPoints[tipKey]
	--print(I_farthestDistance)
		
	--Define the four sets that will be used the most for the rest of the algorithm
	local Facets = {} --{ {k1, k2, k3}, {k2, k3, k4}, etc.} --all representing triangles/facets
	local FacetPointSet = {}   --e.g. FacetPointSet[1] corresponds to points(keys; not actual v3 points) in-front of triangle produced by p1, p2, p3 where pi = SetOfPoints[ki]
	local EdgeSet = {} --want inner sets containing index of triangles sharing the same edge for easy look-up
	--e.g. EdgeSet[k1][k2] = {3, 5} --i.e. k1 and k2 represent setOfPoints key values of the two points that make up the edge. 3 and 5 correspond to facet 3 and 5.
	local InverseEdgeSet = {} --InverseEdgeSet[FacetKeyIndex] = {{k1, k2}, {k2, k3}, {k3, k1}}
	
	--Create initial facets from baseline_p1, baseline_p2, apexPoint, tipPoint for initial simplex
	Facets[1] = {baseline_k1, baseline_k2, apexKey}
	Facets[2] = {baseline_k2, tipKey, apexKey}
	Facets[3] = {apexKey, tipKey, baseline_k1}
	Facets[4] = {tipKey, baseline_k2, baseline_k1}
	--Note to self: from here on out, avoid usage of # and table.insert and table.remove operators (where necessary); if used, comment a star --*
	local function UpdateEdgeSets(NewFacetIndex) --do this after adding a new facet to Facets; will be used to extract horizon edges
		local currentFacet = Facets[NewFacetIndex]
		if InverseEdgeSet[NewFacetIndex] == nil then InverseEdgeSet[NewFacetIndex] = {} end
		for i = 1, 2 do
			for j = i + 1, 3 do
				local k1, k2 = currentFacet[i], currentFacet[j]
				
				if i == 1 and j == 3 then k1, k2 = k2, k1 end
				if EdgeSet[k1] == nil then EdgeSet[k1] = {} end
				if EdgeSet[k2] == nil then EdgeSet[k2] = {} end
				if EdgeSet[k1][k2] == nil then EdgeSet[k1][k2] = {} end 
				if EdgeSet[k2][k1] == nil then EdgeSet[k2][k1] = {} end 
				table.insert(InverseEdgeSet[NewFacetIndex], {k1, k2}) --*
				table.insert(EdgeSet[k1][k2], NewFacetIndex) --*
				table.insert(EdgeSet[k2][k1], NewFacetIndex) --*
			end
		end
	end
	
	
	local initialTetrahedronCentroid = (((baseline_p1 + baseline_p2 + apexPoint)/3) + tipPoint)/2
	local sample_Facet = Facets[1]
	local testInitialNorm = (setOfPoints[sample_Facet[2]] - setOfPoints[sample_Facet[1]]):Cross(setOfPoints[sample_Facet[2]] - setOfPoints[sample_Facet[3]])
	local sample_Centroid = (setOfPoints[sample_Facet[1]] + setOfPoints[sample_Facet[2]] + setOfPoints[sample_Facet[3]])/3
	if testInitialNorm.Unit:Dot((initialTetrahedronCentroid - sample_Centroid).Unit) >= 0 then --prevent normals pointing inwards
		Facets[1] = {apexKey, baseline_k2, baseline_k1}
		Facets[2] = {apexKey, tipKey, baseline_k2}
		Facets[3] = {baseline_k1, tipKey, apexKey}
		Facets[4] = {baseline_k1, baseline_k2, tipKey}
	end
	
	for i = 1, 4 do
		UpdateEdgeSets(i)
		FacetPointSet[i] = {}
	end
	
	
	for currKey, currPoint in pairs(setOfPoints) do --Generate initial facet point sets
		local pointAlreadyAdded = false
		for i = 1, 4 do
			if pointAlreadyAdded == false then
				local currFacet = Facets[i]
				local normVec, posnOnSurf = (setOfPoints[currFacet[2]] - setOfPoints[currFacet[1]]):Cross(setOfPoints[currFacet[2]] - setOfPoints[currFacet[3]]), setOfPoints[currFacet[1]]
				local DFromSurf = DistFromSurface(normVec, posnOnSurf, currPoint)
				if DFromSurf > 0.1 then --beware floating point error, i.e. if points from currFacet are put in it's point set, then no face ever ends up having an empty point set and hence, the stack is always full!
					pointAlreadyAdded = true
					table.insert(FacetPointSet[i], currKey); print(currPoint) --*
				end
			end
		end
	end

	local stack, endOfFacetTableN = {}, 4 --stack contains facet lookup index (e.g. stack[3] = 4 corresponds to Facets[4]
	
	for i = 1, 4 do
		if #FacetPointSet[i] ~= 0 then --if current facet has non-empty point set
			table.insert(stack, i) --*
		end
	end
	
	----Iteration phase
	
	while #stack ~= 0 do --*
		local currentF_index = stack[#stack]
		local currentF = Facets[currentF_index]
		local currentF_farkey, currentF_farpoint = FarthestPointFromSurface(setOfPoints[currentF[1]], setOfPoints[currentF[2]], setOfPoints[currentF[3]], FacetPointSet[currentF_index])
		
		local visibleSet, horizonEdges, visSetHash = {}, {}, {} --contains facet index's of facets that can be seen by currentF_farpoint
		
		--Slightly less efficient method will be used to acquire visibleSet and horizonEdges
		for cKey, cFacet in pairs(Facets) do --build visibleSet
			local p1, p2, p3 = setOfPoints[cFacet[1]], setOfPoints[cFacet[2]], setOfPoints[cFacet[3]]
			local normalVec, posOnSurface = (p2 - p1):Cross(p2 - p3), p1
			local D = DistFromSurface(normalVec, posOnSurface, currentF_farpoint)
			if D > 0 then --beware floating point error
				table.insert(visibleSet, cKey) --*
				visSetHash[cKey] = true --let the hash know that cKey exists in visible set
			end
		end
		
		for i = 1, #visibleSet do --scan all edges in visible set for the boundary edges
			local cKey = visibleSet[i]
			local cEdges = InverseEdgeSet[cKey] -- = {{k1, k2}, {k2, k3}, {k3, k1}}
			for j = 1, 3 do
				local k1, k2 = cEdges[j][1], cEdges[j][2]
				local connectingFacets = EdgeSet[k1][k2] --contains facets containing edge
				for k = 1, #connectingFacets do
					local connectingFacetKey = connectingFacets[k]
					if Facets[connectingFacetKey] ~= nil and visSetHash[connectingFacetKey] == nil then
						--if facet actually exists and it is not in the visible set then the non-visible facet shares an edge with a visible facet
						--i.e. that shared edge must then be on the boundary of the visibleSet
						table.insert(horizonEdges, {k1, k2}) --*
					end
				end
			end
		end
		
		local pointSetHash = {} --hash of visible facets point sets
		
		for i = 1, #horizonEdges do --create new facets
			local cEdge = horizonEdges[i]
			local k1, k2 = cEdge[1], cEdge[2]
			local newFacet = {k2, currentF_farkey, k1}
			endOfFacetTableN = endOfFacetTableN + 1
			Facets[endOfFacetTableN] = newFacet
			UpdateEdgeSets(endOfFacetTableN)
			FacetPointSet[endOfFacetTableN] = {}
			
			for j = 1, #visibleSet do --update facet point sets
				local cKey = visibleSet[j]
				local cPointSet = FacetPointSet[cKey]
				
				for k = 1, #cPointSet do
					local ps_currKey = cPointSet[k]
					if pointSetHash[ps_currKey] == nil then --don't forget to add hash[ps_currkey] = true if in front of new facet
						local normVec, posnOnSurf = (setOfPoints[newFacet[2]] - setOfPoints[newFacet[1]]):Cross(setOfPoints[newFacet[2]] - setOfPoints[newFacet[3]]), setOfPoints[newFacet[1]]
						local DFromSurf = DistFromSurface(normVec, posnOnSurf, setOfPoints[ps_currKey])
						
						if DFromSurf > 0.1 then
							pointSetHash[ps_currKey] = true
							table.insert(FacetPointSet[endOfFacetTableN], ps_currKey) --*
						end
					end
				end
			end
			
			if #FacetPointSet[endOfFacetTableN] ~= 0 then --if current new facet has non-empty point set
				table.insert(stack, endOfFacetTableN) --*
			end
		end

		for i = 1, #visibleSet do --remove faces in visible set from facets and from stack
			local removeFacetKey = visibleSet[i]
			Facets[removeFacetKey] = nil
			
			for j = 1, #stack do
				local theCurrentStack = stack[j]
				if theCurrentStack == removeFacetKey then
					table.remove(stack, j)
				end
			end
		end
	end

	return Facets
end
