#' Interactively tune model parameters with a Shiny gadget
#'
#' Launches a small Shiny gadget with sliders in the left sidebar for each
#' numeric element of `pars` and a live-updating plot of
#' `plotAgeLikelihood(pars, age_at_length)` in the main panel.
#'
#' @param pars A named list of scalar numeric parameters. Typical entries are
#'   `k`, `L_inf`, `d`, `m`, `spawning_mu`, `spawning_kappa`, `annuli_date`,
#'   and `annuli_min_age`.
#' @param age_at_length A data frame of observed age-at-length data as used by
#'   `plotAgeLikelihood()` (see `Cod_CS_age_at_length`).
#'
#' @return Invisibly returns the final parameters after the gadget is closed.
#'   While the gadget is running it displays the plot interactively.
#'
#' @examples
#' # Not run: requires interactive session
#' # data("Cod_CS_age_at_length")
#' # data("Cod_CS_pars")
#' # tune_pars(Cod_CS_pars, Cod_CS_age_at_length)
#'
#' @export
tune_pars <- function(pars, age_at_length) {

	stopifnot(is.list(pars))

	# Keep only scalar numerics
	numeric_pars <- pars[vapply(pars, function(x) is.numeric(x) && length(x) == 1, logical(1))]
	if (length(numeric_pars) == 0) stop("No scalar numeric parameters found in `pars`.")

	# Helper to pick slider ranges based on name and value
	make_slider <- function(name, value) {
		value_num <- as.numeric(value)
		# Heuristics for common parameters
		if (grepl("^(spawning_mu|annuli_date)$", name)) {
			minv <- 0; maxv <- 1; step <- 0.001
		} else if (grepl("^(spawning_kappa)$", name)) {
			minv <- 0; maxv <- max(5, value_num * 3); step <- 0.1
		} else if (grepl("^(annuli_min_age)$", name)) {
			minv <- 0; maxv <- max(2, value_num * 3); step <- 0.01
		} else if (grepl("^L_inf$", name)) {
			minv <- max(0, value_num * 0.25);
			maxv <- max(minv + 1, value_num * 2)
			step <- max(0.1, round((maxv - minv) / 200, 2))
		} else {
			# Generic non-negative numeric
			minv <- 0
			maxv <- if (value_num > 0) value_num * 3 else 1
			step <- max(10^floor(log10(max(1e-6, maxv))/ -2), 1e-4)
		}
		shiny::sliderInput(inputId = name, label = name,
						 min = minv, max = maxv, value = value_num, step = step)
	}

	ui <- shiny::fluidPage(
		shiny::titlePanel("Tune parameters"),
		shiny::sidebarLayout(
			shiny::sidebarPanel(
				shiny::tagList(
					lapply(names(numeric_pars), function(nm) make_slider(nm, numeric_pars[[nm]]))
				)
			),
			shiny::mainPanel(
				shiny::plotOutput("lik_plot", height = "600px")
			)
		)
	)

	server <- function(input, output, session) {
		current_pars <- shiny::reactive({
			updated <- pars
			for (nm in names(numeric_pars)) {
				updated[[nm]] <- input[[nm]]
			}
			updated
		})

		output$lik_plot <- shiny::renderPlot({
			plotAgeLikelihood(current_pars(), age_at_length)
		})

		# Return updated parameters when session ends
		shiny::onSessionEnded(function() {
			invisible(current_pars())
		})
	}

	# Run as gadget so it appears as a small window in RStudio, but fall back to app
	app <- shiny::shinyApp(ui, server)
	if (requireNamespace("miniUI", quietly = TRUE)) {
		shiny::runGadget(app, viewer = shiny::dialogViewer("Tune parameters", width = 1100, height = 800))
	} else {
		shiny::runApp(app)
	}
}


