### FUNCTION TO IMPORT DATA
### DELETE ENTRIES BEYOND THE  ZERO 'BASELINE'
function import_data(fileno, MTT)
	extn = string("tendon_data/mtt",MTT,"/converted_data/stress_strain",fileno,".txt")
	global data = readdlm(extn, ',')
	# delete non-unique entries (w.r.t strain)
	global row_count = 2;

	while (row_count <= size(data)[1])
		if data[row_count,1] == data[row_count-1,1]
			global data = data[1:end .!=row_count, :]
		else
			global row_count += 1
		end
	end

	# convert to MPa and stretch
	data[:,2] = data[:,2]/10^6
	data[:,1] = data[:,1] .+ 1


	# compute baseline and delete entries beyond
	baseline = data[end,2]
	max_stress = maximum(data[:,2])
	tol = 0.1*max_stress

	I = findall(I->I==max_stress, data[:,2])[end]

	cutoff = length(data[:,2])
	for i=I:length(data[:,2])
		if abs(data[i,2]-baseline) < tol
			cutoff = i;
			break;
		end
	end

	# define new arrays for the non-zero data values
	nonzero_strain = data[1:cutoff,1]
	nonzero_stress = data[1:cutoff,2]
	nonzero_data = [nonzero_strain nonzero_stress]

	return nonzero_data

end

# import the data as above, but return the baseline
function import_baseline(fileno, MTT)
	extn = string("tendon_data/mtt",MTT,"/converted_data/stress_strain",fileno,".txt")
	global data = readdlm(extn, ',')
	# delete non-unique entries (w.r.t strain)
	global row_count = 2;

	while (row_count <= size(data)[1])
		if data[row_count,1] == data[row_count-1,1]
			global data = data[1:end .!=row_count, :]
		else
			global row_count += 1
		end
	end

	data[:,2] = data[:,2]/10^6
	data[:,1] = data[:,1] .+ 1


	# compute baseline and delete entries beyond
	baseline = data[end,2]
	return baseline

end




### DEFINE A TRIANGULAR DISTRIBUTION
function triangular_dist(x, a, b, c)
	# account for the delta function case..
	# hacky
	if (a == b)
		b += 0.0001
	end

	l = length(x)
	f = ones(l)
	for i = 1:l
		if (a <= x[i]) && (x[i] < c)
			f[i] = 2*(x[i]-a)/((b-a)*(c-a))
		elseif (c <= x[i]) && (x[i] < b)
			f[i] = 2*(b-x[i])/((b-a)*(b-c))
		else
			f[i] = 0
		end
	end
	return f
end

### BIMODAL TRIANGULAR DISTRIBUTION DEFINITION
function bi_modal_triangular_dist(x, a1, b1, c1, a2, b2, c2, M2)
	return (triangular_dist(x,a1,b1,c1) + M2*triangular_dist(x,a2,b2,c2))/(1+M2)
end

### ELASTIC STRESS IN A FIBRIL - FOR USE WITH THE ELASTIC MODEL ONLY
function fibril_stress_elastic(x, lambda_C, E)
	l = length(lambda_C)
	f = ones(l)
	for i = 1:l
		if x <= lambda_C[i]
			f[i] = 0
		elseif x > lambda_C[i]
			f[i] = E*(x/lambda_C[i] - 1)
		end
	end
	return f
end

### ELASTOPLASTIC STRESS IN FIBRIL
function fibril_stress_elastoplastic(x, lambda_C, lambda_Y, lambda_R, E, k)
	 l = length(lambda_C)
	 f = ones(l)
	 for i = 1:l
	     if (x > lambda_C[i]) && (x <= lambda_C[i]*lambda_Y[i])
	     	f[i] = E*(x/lambda_C[i] - 1)
		elseif (x > lambda_C[i]*lambda_Y[i]) && (x <= lambda_C[i]*lambda_R[i])
	     	f[i] = E*(1-k)*(lambda_Y[i]-1) + E*k*(x/lambda_C[i] - 1)
		else
			f[i] = 0
	     end

	 end

	 return f
end


### COMPUTE ELASTIC TENDON STRESS USING INTEGRAL METHOD
function elastic_tendon_stress(x, a, b, c, E, phi)
	l = length(x)
	stress = zeros(l)
	for i = 1:l
		integrand(X) = fibril_stress_elastic(x[i],X,E).*triangular_dist(X, a, b, c);
		stress[i] = hquadrature(integrand, 1, x[i])[1][1]
	end

	return phi*stress
end


### ANALYTIC EXPRESSION FOR ELASTIC TENDON STRESS
function elastic_tendon_stress_analytic(x, a, b, c, E, phi)

	l = length(x)
	stress = zeros(l)

	for i = 1:l
		if  (x[i] >= a) && (x[i] < c)
			stress[i] += (-a^2)/((a-b)*(a-c)) + x[i]*(2*a*log(a))/((a-b)*(a-c)) +
			(x[i]^2)/((a-b)*(a-c)) + (x[i]*log(x[i]))*(-2*a)/((a-b)*(a-c));
		elseif (x[i] >= c) && (x[i] < b)
			stress[i] += (-b*c)/((a-b)*(b-c))-a/(a-b) +
			x[i]*((2*a*log(a))/((a-b)*(a-c)) + (2*c*log(c))/((a-c)*(b-c))) +
			(x[i]^2)/((a-b)*(b-c)) + (x[i]*log(x[i]))*(-2*b)/((a-b)*(b-c));
		elseif x[i] > b
			stress[i] += -1 +
			x[i]*((2*a*log(a))/((a-b)*(a-c)) + (2*c*log(c))/((a-c)*(b-c)) +
			(-2*b*log(b))/((a-b)*(b-c))) ;
		else
			stress[i] += 0;
		end
	end

	return phi*E*stress;

end

# introduce penalization function for the analytic stress
# to be used in the asymmetric case
function elastic_tendon_stress_analytic_penalised(x, a, b, c, E, phi, lambdaY)
	if (b > lambdaY)
		penalise = 10000
	else
		penalise = 0
	end

	return elastic_tendon_stress_analytic(x, a, b, c, E, phi) .+ penalise
end


### COMPUTE TENDON STRESS USING 3D INTEGRAL
### ALL STRETCH DISTRIBUTIONS ARE ASYMMETRIC
### AND TRIANGULAR.
function tendon_stress_symmetric_unimodal(stretch, aY, bY, aR, bR, k, aC, bC, E, phi, TOL)
	#println("pars = ", [ aY, bY, aR, bR])
	l = length(stretch)
	stress=zeros(l)

	cC = (aC + bC)/2.0
	cY = (aY + bY)/2.0
	cR = (aR + bR)/2.0

	for i = 1:l
		integrand3(x) = phi*fibril_stress_elastoplastic(stretch[i], x[1], x[2], x[3], E, k).*
		triangular_dist(x[1], aC, bC, cC).*
		triangular_dist(x[2], aY, bY, cY).*
		triangular_dist(x[3], aR, bR, cR);
		stress[i] = hcubature(integrand3, [1, aY, aR],
		[stretch[i], bY, bR]; rtol=TOL, initdiv = 4, maxevals=10000)[1][1]
	end
	#println("Stress computed")
	return stress
end


function tendon_stress_symmetric_unimodal_split(stretch, aY, bY, aR, bR, k, aC, bC, E, phi, TOL)
	#println("pars = ", [ aY, bY, aR, bR])
	l = length(stretch)
	stress=zeros(l)

	cC = (aC + bC)/2.0
	cY = (aY + bY)/2.0
	cR = (aR + bR)/2.0

	# split the integrand up into multiple sections - convinient because
	# the distribution functions are piecewise linear.
	for i = 1:l
		integrand3(x) = phi*fibril_stress_elastoplastic(stretch[i], x[1], x[2], x[3], E, k).*
		triangular_dist(x[1], aC, bC, cC).*
		triangular_dist(x[2], aY, bY, cY).*
		triangular_dist(x[3], aR, bR, cR);

		if (stretch[i]) >= 1 && (stretch[i] < cC)
			stress[i] += hcubature(integrand3, [1, aY, aR],
			[stretch[i], cY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, aR],
			[stretch[i], bY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, aY, cR],
			[stretch[i], cY, bR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, cR],
			[stretch[i], bY, bR]; rtol=TOL)[1][1]
		end

		if (stretch[i]) >= cC && (stretch[i] < bC)
			stress[i] += hcubature(integrand3, [1, aY, aR],
			[cC, cY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, aR],
			[cC, bY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, aY, cR],
			[cC, cY, bR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, cR],
			[cC, bY, bR]; rtol=TOL)[1][1]

			stress[i] += hcubature(integrand3, [cC, aY, aR],
			[stretch[i], cY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, cY, aR],
			[stretch[i], bY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, aY, cR],
			[stretch[i], cY, bR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, cY, cR],
			[stretch[i], bY, bR]; rtol=TOL)[1][1]
		end

		if stretch[i] >= bC
			stress[i] += hcubature(integrand3, [1, aY, aR],
			[cC, cY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, aR],
			[cC, bY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, aY, cR],
			[cC, cY, bR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [1, cY, cR],
			[cC, bY, bR]; rtol=TOL)[1][1]

			stress[i] += hcubature(integrand3, [cC, aY, aR],
			[bC, cY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, cY, aR],
			[bC, bY, cR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, aY, cR],
			[bC, cY, bR]; rtol=TOL)[1][1]
			stress[i] += hcubature(integrand3, [cC, cY, cR],
			[bC, bY, bR]; rtol=TOL)[1][1]
		end

	end
	return stress
end




### SAME AS ABOVE BUT WITH SYMMETRIC DISTRIBUTIONS
function tendon_stress_symmetric_bimodal(stretch, aY, bY, aR1, bR1, aR2, bR2, M, k, aC, bC, E, phi, TOL)

	l = length(stretch)
	stress=zeros(l)

	cC = (aC + bC)/2.0
	cY = (aY + bY)/2.0
	cR1 = (aR1 + bR1)/2.0
	cR2 = (aR2 + bR2)/2.0

	for i =1:l
		integrand3(x) = phi*fibril_stress_elastoplastic(stretch[i], x[1], x[2], x[3], E, k).*
		triangular_dist(x[1], aC, bC, cC).*
		triangular_dist(x[2], aY, bY, cY).*
		bi_modal_triangular_dist(x[3], aR1, bR1, cR1, aR2, bR2, cR2, M);
		stress[i] = hcubature(integrand3, [1, aY, min(aR1,aR2)],
		[stretch[i], bY, max(bR1,bR2)]; rtol=TOL, initdiv=4, maxevals=10000)[1][1]
	end

	return stress
end

function tendon_stress_symmetric_unimodal_penalised(stretch, aY, bY, aR, bR, k, aC, bC, E, phi, TOL)
	penalty = 0
	if (bY <= aY)||(bR <= aR)
		penalty = 10000
	end
	return tendon_stress_symmetric_unimodal(stretch, aY, bY, aR, bR, k, aC, bC, E, phi, TOL) .+ penalty
end


function tendon_stress_symmetric_bimodal_penalised(stretch, aY, bY, aR1, bR1, aR2, bR2, M, k, aC, bC, E, phi, TOL)
	penalty = 0
	if (bY <= aY)||(bR1 <= aR1)||(bR2 <= aR2)
		penalty = 10000
	end
	return tendon_stress_symmetric_bimodal(stretch, aY, bY, aR1, bR1, aR2, bR2, M, k, aC, bC, E, phi, TOL) .+ penalty
end

### FITTING FUNCTION FOR THE ELASTIC CASE
### USE A GLOBAL SYMMETRY VARIABLE TO SELECT
### BETWEEN SYMMETRIC AND ASYMMETRIC TRIANGULAR
### DISTRIBUTIONS
function fit_elastic_data(fileno, yield_stretch, phi, MTT, io)
	# first get relevant data
	data = import_data(fileno, MTT)
	# isolate the elastic data
	cutoff = 4
	for i = 1:length(data)
		if data[i,1] > yield_stretch
			if i >= 4
				cutoff = i;
			end
			break;
		end
	end
	elastic_data = data[1:cutoff,:]


	xdata = elastic_data[:,1]
	ydata = elastic_data[:,2]

	l = length(xdata)
	residuals = zeros(l)

	p0 = [(yield_stretch-1)*(2/3)+1, 500]
	lb = [1.0, 0.0]
	ub = [yield_stretch, 5000]

	# ef = elastic_fit
	ef = so.curve_fit((x, b, E)->elastic_tendon_stress_analytic_penalised(x, 1.0, b, (b+1.0)/2.0, E, phi, yield_stretch),
	xdata, ydata, p0, bounds = (lb, ub))[1]

	# calculate residuals - just use analytic fit here
	residuals = ydata - elastic_tendon_stress_analytic(xdata, 1.0, ef[1], (1.0+ef[1])/2.0, ef[2], phi)

	# compute RMSE
	RMSE = sqrt(mean(residuals.^2))

	# print results of fitting to terminal
	println(string("Elastic fitting parameters: b = ", ef[1] ,", E = ", ef[2]))
	# write the results to a file
	write(io, @sprintf "%s %d %.2f %.3f %.3f %.3f %.3f \n" MTT fileno phi yield_stretch ef[1] ef[2] RMSE)

	return ef
end

# function to check if the data should be excluded
function exclude_data(MTT, file_no, yield_stretch, io)

	# first get relevant data
	data = import_data(file_no, MTT)
	cutoff = 0
	# isolate the elastic data
	for i = 1:length(data)
		if data[i,1] > yield_stretch
			cutoff = i
			break
		end
	end
	elastic_data = data[1:cutoff,:]
	plastic_data = data[cutoff + 1:end,:]

	baseline = import_baseline(file_no, MTT)

	if (baseline > 0.5*maximum(data[:,2]))
		println("Non-zero final stress: exclude from fitting")
		println(io, @sprintf "Exclude %s %d , non-zero final stress" MTT file_no)
		return true
	end

	if (cutoff < 4)
		println("Fewer than 4 elastic points: exclude from fitting")
		println(io, @sprintf "Exclude %s %d , less than 4 elastic points" MTT file_no)
		return true
	end

	if (length(plastic_data[:,1]) < 7)
		println("Fewer than 7 plastic points: exclude from fitting")
		println(io, @sprintf "Exclude %s %d , less than 7 plastic points" MTT file_no)
		return true
	end


	return false
end

# function for determining a random initial guess
function random_p0(peak_stretch, max_stretch, yield_stretch)
	p0_1 = yield_stretch
	p0_2 = rand(Uniform(p0_1, max_stretch))
	p0_3 = rand(Uniform(p0_1, max_stretch))
	p0_4 = rand(Uniform(p0_3, max_stretch))
	p0_5 = 0

	return [p0_1, p0_2, p0_3, p0_4, p0_5]

end

function random_bimodal_p0(peak_stretch, max_stretch, yield_stretch)


	p0_1 = yield_stretch	# aY
	p0_2 = rand(Uniform(p0_1, max_stretch)) # bY
	p0_3 = rand(Uniform(p0_1, max_stretch)) # aR1
	p0_4 = rand(Uniform(p0_3, max_stretch)) # bR1
	p0_5 = rand(Uniform(p0_1, max_stretch)) # aR2
	p0_6 = rand(Uniform(p0_5, max_stretch)) # bR2
	p0_7 = 0.5	# M
	p0_8 = 0.0	# k

	return [p0_1, p0_2, p0_3, p0_4, p0_5, p0_6, p0_7, p0_8]

end


### FITTING FUNCTION FOR THE FULL TENDON STRESS
function fit_plastic_data_unimodal(fileno, yield_stretch, phi, MTT, bC, E, lb, ub, p0, TOL)
	# first get relevant data
	data = import_data(fileno, MTT)

	xdata = data[:,1]
	ydata = data[:,2]

	l = length(xdata)
	residuals = zeros(l)

	## maximum value of strain in the data
	max_x = maximum(xdata)
	max_y = maximum(ydata)

	pf = so.curve_fit((x, p1, p2, p3, p4, p5)->
	tendon_stress_symmetric_unimodal_penalised(x, p1, p2, p3, p4, p5, 1.0, bC, E, phi, TOL),
	xdata, ydata, p0, bounds = (lb, ub))[1]

	# compute residuals
	residuals = ydata - tendon_stress_symmetric_unimodal(xdata, pf[1], pf[2], pf[3], pf[4], pf[5], 1.0, bC, E, phi, TOL)

	# compute RMSE
	RMSE = sqrt(mean(residuals.^2))

	return[pf, RMSE]

end


### FITTING FUNCTION FOR THE BIMODAL TENDON STRESS
function fit_plastic_data_bimodal(fileno, yield_stretch, phi, MTT, bC, E, lb, ub, p0, TOL)
	# first get relevant data
	data = import_data(fileno, MTT)

	xdata = data[:,1]
	ydata = data[:,2]

	l = length(xdata)
	residuals = zeros(l)

	## maximum value of strain in the data
	max_x = maximum(xdata)
	max_y = maximum(ydata)

	pf = so.curve_fit((x, p1, p2, p3, p4, p5, p6, p7, p8)->
	tendon_stress_symmetric_bimodal_penalised(x, p1, p2, p3, p4, p5, p6, p7, p8, 1.0, bC, E, phi, TOL),
	xdata, ydata, p0, bounds = (lb, ub))[1]

	# compute residuals
	residuals = ydata - tendon_stress_symmetric_bimodal(xdata, pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], 1.0, bC, E, phi, TOL)

	# compute RMSE
	RMSE = sqrt(mean(residuals.^2))

	return[pf, RMSE]

end

# function to check if the data should be excluded
function exclude_data_bimodal(MTT, file_no, yield_stretch, io)

	# first get relevant data
	data = import_data(file_no, MTT)
	cutoff = 0
	# isolate the elastic data
	for i = 1:length(data)
		if data[i,1] > yield_stretch
			cutoff = i
			break
		end
	end
	plastic_data = data[cutoff + 1:end,:]

	if (length(plastic_data[:,1]) < 10)
		println("Fewer than 10 plastic points: cannot fit using bimodal rupture distribution")
		println(io, @sprintf "Exclude %s %d from bimodal fits, less than 10 plastic points" MTT file_no)
		return true
	end


	return false
end
