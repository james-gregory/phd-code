using Plots; pyplot()
using QuadGK
using HCubature
using LsqFit
using DelimitedFiles
using PyCall
using Printf
using DataFrames
using CSV
using Statistics
using Random
using Distributions
using BenchmarkTools

# include the function file
include("functions.jl")

# import a python library for curve fitting
so = pyimport("scipy.optimize")

# For convinience, wrap everything up into a function and specify start and end
# points d1 and d2. These are the data indices that we start and end with.
function fit_data(d1, d2)
	# PUT METADATA HERE
	# import all metadata
	df = CSV.read("GOH.csv", DataFrame)

	# get number of rows
	n_row = nrow(df)

	# prepare files for output
	io_elastic = open(string("output/elastic_pars_", d1, "-", d2, ".txt"), "a")
	io_plastic = open(string("output/plastic_pars_", d1, "-", d2, ".txt"), "a")
	io_bimodal = open(string("output/bimodal_pars_", d1, "-", d2, ".txt"), "a")
	io_initial = open(string("output/initial_guess_", d1, "-", d2, ".txt"), "a")
	io_ER = open(string("output/elastic_rupture_pars_", d1, "-", d2, ".txt"), "a")
	io_exclude = open(string("output/exclude_", d1, "-", d2, ".txt"), "a")

	# loop over data
	for i = d1:d2

		#extract metadata
		if df[i, "MTT"] < 10
			mtt = string("0", df[i, "MTT"])
		else
			mtt = string(df[i, "MTT"])
		end

		file_no = Int(df[i, "stress_strain"])
		yield_stretch = df[i, "Yield stretch"]
		phi = df[i, "phi"]

		steps = df[i, "Steps"]

		# output the test number to keep track
		@printf "Test: %s - %d \n" mtt file_no

		# check if the data should be excluded or not
		if(exclude_data(mtt, file_no, yield_stretch, io_exclude))
			continue
		end

		# do the elastic fitting
		ef = fit_elastic_data(file_no,yield_stretch,phi, mtt, io_elastic)

		# import the data
		plot_data = import_data(file_no, mtt)
		max_stretch = maximum(plot_data[:,1])
		peak_stretch = plot_data[findmax(plot_data[:,2])[2],1]

		# set boundaries for the fitting
		lb = [1, 1, 1, 1, 0.0]
		ub = [peak_stretch, max_stretch, max_stretch, max_stretch, 1.0]

		# perform multiple fits with random initial guesses
		# initialise a p0
		p0_best = [0, 0, 0, 0, 0]

		# initialise an error
		err = 1000
		bimodal_err = 1000
		no_attempts = 50
		max_time = 60
		wait_time = 5
		global current_fit = 0
		global best_fit = 0
		global current_bimodal_fit = 0
		global best_bimodal_fit = 0
		j = 1
		while (j <= no_attempts)

			@printf "Attempt %d/%d. Current error = %.3f. \n" j no_attempts err
			# extract some random p0 (low tol)
			p0_current = random_p0(peak_stretch, max_stretch, yield_stretch)

			try current_fit = fit_plastic_data_unimodal(file_no, yield_stretch, phi, mtt, ef[1], ef[2], lb, ub, p0_current, 0.001)
				j += 1
			catch
				println("Parameters not found - will try again.")
				continue
			end

			# if the error is lower with this initial guess
			if current_fit[2] < err
				# set new lowest error
				err = current_fit[2]
				# set new best initial guess
				p0_best = p0_current
				# set the new best fit (bf)
				best_fit = current_fit
			end
		end

		bf = best_fit
		# output initial guess to file
		write(io_initial, @sprintf "%.3f %.3f %.3f %.3f %.3f \n" p0_best[1] p0_best[2] p0_best[3] p0_best[4] p0_best[5])
		# print best fit parameters
		println("Best fit parameters: ", bf[1])
		# output the best fit parameters to a file
		write(io_plastic, @sprintf "%s %d %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n" mtt file_no phi yield_stretch bf[1][1] bf[1][2] bf[1][3] bf[1][4] bf[1][5] bf[2])


		### CHECK IF THERE ARE STEPS, AND IF SO REDO WITH BIMODAL RUPTURE DIST
		### also check there's enough points to do the bimodal fit
		if (steps == 1) && !exclude_data_bimodal(mtt, file_no, yield_stretch, io_exclude)
			bimodal = true
		else
			bimodal = false
		end

		if bimodal
			println("Performing bimodal fit.")
			# then do fit using previous best fit parameters as initial guess
			# aY / bY / aR1 / bR1 / aR2 / bR2 / M / k
			lbb = [1, 1, 1, 1, 1, 1, 0.0, 0.0]
			ubb = [peak_stretch, max_stretch, max_stretch, max_stretch, max_stretch, max_stretch, 1.0, 1.0]

			k = 1
			while (k <= no_attempts)

				@printf "Bimodal attempt %d/%d. Current error = %.3f. \n" k no_attempts bimodal_err
				# extract some random p0 (low tol)
				p0_current_bimodal = random_bimodal_p0(peak_stretch, max_stretch, yield_stretch)

				try current_bimodal_fit = fit_plastic_data_bimodal(file_no, yield_stretch, phi, mtt, ef[1], ef[2], lbb, ubb, p0_current_bimodal, 0.001)
					k += 1
				catch
					println("Bimodal parameters not found - will try again.")
					continue
				end

				# if the error is lower with this initial guess
				if current_bimodal_fit[2] < bimodal_err
					# set new lowest error
					bimodal_err = current_bimodal_fit[2]
					# set new best initial guess
					p0_best_bimodal = p0_current_bimodal
					# set the new best fit (bf)
					best_bimodal_fit = current_bimodal_fit
				end
			end

			bimodal_fit = best_bimodal_fit

			println("Bimodal parameters = ", bimodal_fit[1])
			@printf "Error changed from %.2f to %.2f. \n" bf[2] bimodal_fit[2]
			# output to file
			write(io_bimodal, @sprintf "%s %d %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n" mtt file_no phi yield_stretch bimodal_fit[1][1] bimodal_fit[1][2] bimodal_fit[1][3] bimodal_fit[1][4] bimodal_fit[1][5] bimodal_fit[1][6] bimodal_fit[1][7] bimodal_fit[1][8] bf[2])
		end

		# initialise vectors for plotting
		X = LinRange(1,1.0 + 1.1*(max_stretch-1.0), 250)
		Fe = zeros(length(X))
		#Fe = elastic_tendon_stress_analytic(X, 1, ef[1], (1.0+ef[1])/2.0, ef[2], phi)
		Fp = tendon_stress_symmetric_unimodal(X, bf[1][1], bf[1][2], bf[1][3], bf[1][4], bf[1][5], 1.0, ef[1], ef[2], phi, 0.01)



		# plot both fits (and any other functions is necessary)
		scatter(plot_data[:,1],plot_data[:,2])
		# plot ER fit
		#plot!(X,Fe)
		plot!(X,Fp)

		# also plot the bimodal case if it exists
		if bimodal
			Fbimodal = zeros(length(X))
			Fbimodal = tendon_stress_symmetric_bimodal(X, bimodal_fit[1][1], bimodal_fit[1][2], bimodal_fit[1][3], bimodal_fit[1][4], bimodal_fit[1][5], bimodal_fit[1][6], bimodal_fit[1][7], bimodal_fit[1][8], 1.0, ef[1], ef[2], phi, 0.01)
			plot!(X,Fbimodal)
		end

		savefig(string("Figs/fitting_fig_",mtt,"-",file_no,".pdf"))

		flush(io_elastic)
		flush(io_plastic)
		flush(io_bimodal)
		flush(io_ER)
		flush(io_exclude)
		flush(io_initial)


	end

	close(io_elastic)
	close(io_plastic)
	close(io_bimodal)
	close(io_ER)
	close(io_exclude)
	close(io_initial)


end
