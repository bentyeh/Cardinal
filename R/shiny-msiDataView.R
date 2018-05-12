
#### View MSI data in a shinydashboard ####
## ---------------------------------------

msiDataViewUI <- function(id, label = "Data View") {
	ns <- NS(id)
	verticalLayout(
		actionButton(ns("close"), "Close",
			style="float:right; padding:0px 2px 0px 2px; font-size:90%;"),
		uiOutput(ns("msiDataView"))
	)
}

msiDataView <- function(input, output, session, dataset) {

	#### Session variables
	
	ns <- session$ns

	data <- reactive({
		tryCatch(get(dataset, envir=globalenv()),
			error=function(e) NULL)
	})

	sv <- list(
		mz = syncVal(mz(data())[1], function(mz) {
			validate(need(mz, "Invalid m/z value"))
			mz(data())[features(data(), mz=mz)]
		}),
		mz_tol = syncVal(0.001),
		xy = syncVal(as.numeric(coord(data())[1,c(1,2)])),
		xy_names = syncVal(names(coord(data()))[c(1,2)]),
		coord_names = syncVal(coordLabels(data())),
		ionimage_xylim = syncVal(
			c(range(coord(data())[,1]),
				range(coord(data())[,2]))
		),
		spectrum_massrange = syncVal(range(mz(data()))),
		ionimage_intensity_range = syncVal(NULL),
		spectrum_intensity_range = syncVal(NULL),
		ionimage_contrast = syncVal("none"),
		ionimage_smoothing = syncVal("none"),
		ionimage_colorscale = syncVal("rainbow"),
		ionimage_function = syncVal("mean"),
		plot_layout = syncVal(50), # 50%, 50%
		ionimage_height = syncVal(400), # 400px
		spectrum_height = syncVal(400), # 400px
		closed = reactiveVal(FALSE)
	)

	sv[["mz_range"]] <- reactive({
		range(mz(data()))
	})

	sv[["xy_range"]] <- reactive({
		c(range(coord(data())[,sv$xy_names()[1]]),
			range(coord(data())[,sv$xy_names()[2]]))
	})

	sv[["subset_all"]] <- reactive({
		xy <- sv$xy_names()
		coord <- names(coord(data()))
		sub <- !coord %in% xy
		if ( any(sub) ) {
			coord <- coord(data())[coord[sub]]
			sub <- expand.grid(lapply(coord, unique))
		} else {
			NULL
		}
	})

	sv[["subset_names"]] <- reactive({
		if ( !is.null(sv$subset_all()) ) {
			apply(sv$subset_all(), 1, function(a) {
				paste0(paste0(names(a), "=='", a, "'"), collapse=" & ")
			})
		} else {
			""
		}
	})

	sv[["subset"]] <- syncVal(sv$subset_names()[1])

	sv[["data"]] <- data

	sv[["name"]] <- dataset

	### Helper functions

	# Get logical subset
	get_subset_logical <- reactive({
		if ( sv$subset() == "" ) {
			TRUE
		} else {
			eval(parse(text=sv$subset()), as.list(coord(data())))
		}
	})

	# Get position subset
	get_subset_pos <- reactive({
		if ( sv$subset() == "" ) {
			NULL
		} else {
			idx <- which(sv$subset() == sv$subset_names())
			sv$subset_all()[idx,,drop=FALSE]
		}
	})

	# close view
	observeEvent(input$close, {
		if ( input$close > 0 )
			sv$closed(TRUE)
	})

	#### Plot output

	output$msiViewer <- renderUI({
		validate(need(sv$plot_layout(), "Invalid plot layout"))
		validate(need(sv$ionimage_height(), "Invalid plot height"))
		validate(need(sv$spectrum_height(), "Invalid plot height"))
		splitLayout(
			cellWidths=paste0(c(sv$plot_layout(), 100 - sv$plot_layout()), "%"),
			plotOutput(ns("ionimage"),
				click=clickOpts(id=ns("ionimage_click")),
				dblclick=dblclickOpts(id=ns("ionimage_dblclick")),
				brush=brushOpts(id=ns("ionimage_brush"),
					direction="xy", resetOnNew=TRUE),
				height=paste0(sv$ionimage_height(), "px")),
			plotOutput(ns("spectrum"),
				click=clickOpts(id=ns("spectrum_click")),
				dblclick=dblclickOpts(id=ns("spectrum_dblclick")),
				brush=brushOpts(id=ns("spectrum_brush"),
					direction="x", resetOnNew=TRUE),
				height=paste0(sv$spectrum_height(), "px"))
		)
	})

	# outputOptions(output, ns("ionimage"), priority=-1)
	# outputOptions(output, ns("spectrum"), priority=-1)

	plot_null <- function() {
		par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0),
			cex.axis=1, cex.lab=1)
		plot(0, 0, type='n', xlab="", ylab="",
			xaxt='n', yaxt='n')
		text(0, 0, "Nothing to plot.")
	}

	plot_spectrum <- function() {
		validate(
			need(sv$xy(), "Invalid x/y position"),
			need(sv$xy_names(), "Invalid x/y names"),
			need(sv$spectrum_massrange(), "Invalid mass range"),
			need(!anyNA(sv$spectrum_intensity_range()), "Invalid intensity range")
		)
		par(mar=c(5,3,2,1))
		pos <- list(sv$xy()[1], sv$xy()[2])
		names(pos) <- sv$xy_names()
		pos <- c(pos, get_subset_pos())
		fm <- paste0("~", "mz")
		tryCatch({
			if ( is.null(sv$spectrum_intensity_range()) ) {
				plot(sv$data(),
					as.formula(fm),
					coord=pos,
					xlim=sv$spectrum_massrange())
			} else {
				plot(sv$data(),
					as.formula(fm),
					coord=pos,
					xlim=sv$spectrum_massrange(),
					ylim=sv$spectrum_intensity_range())
			}
		}, warning=function(e) plot_null())
	}

	plot_ionimage <- function() {
		validate(
			need(sv$mz(), "Invalid m/z value"),
			need(sv$mz_tol(), "Invalid m/z tolerance"),
			need(sv$xy_names(), "Invalid x/y names"),
			need(sv$ionimage_xylim(), "Invalid x/y limits"),
			need(!anyNA(sv$ionimage_intensity_range()), "Invalid intensity range")
		)
		colors <- color.map(sv$ionimage_colorscale())
		par(mar=c(5,3,2,1))
		fm <- paste0("~", paste0(sv$xy_names(), collapse="*"))
		tryCatch({
			if ( is.null(sv$ionimage_intensity_range()) ) {
				image(sv$data(),
					as.formula(fm),
					mz=sv$mz(),
					plusminus=sv$mz_tol(),
					contrast.enhance=sv$ionimage_contrast(),
					smooth.image=sv$ionimage_smoothing(),
					xlim=sv$ionimage_xylim()[c(1,2)],
					ylim=sv$ionimage_xylim()[c(3,4)],
					fun=match.fun(sv$ionimage_function()),
					col.regions=colors,
					subset=get_subset_logical())
			} else {
				image(sv$data(),
					as.formula(fm),
					mz=sv$mz(),
					plusminus=sv$mz_tol(),
					contrast.enhance=sv$ionimage_contrast(),
					smooth.image=sv$ionimage_smoothing(),
					xlim=sv$ionimage_xylim()[c(1,2)],
					ylim=sv$ionimage_xylim()[c(3,4)],
					zlim=sv$ionimage_intensity_range(),
					fun=match.fun(sv$ionimage_function()),
					col.regions=colors,
					subset=get_subset_logical())
			}
		}, warning=function(e) plot_null())
	}

	plot_mz_marker <- function() {
		mz <- c(sv$mz() - sv$mz_tol(), sv$mz() + sv$mz_tol())
		rect(mz[1], par("usr")[3], mz[2], par("usr")[4],
			col=rgb(1, 0, 0, 0, alpha=0.5), border=NA)
		abline(v=sv$mz(), lty=2, lwd=2, col="blue")
	}

	plot_pos_marker <- function() {
		points(sv$xy()[1], sv$xy()[2], pch=4, lwd=4, cex=2, col="black")
		points(sv$xy()[1], sv$xy()[2], pch=4, lwd=2, cex=2, col="white")
	}

	#### Plot reactivity

	output$spectrum <- renderPlot({
		plot_spectrum()
		plot_mz_marker()
	}, bg="transparent")

	output$ionimage <- renderPlot({
		plot_ionimage()
		plot_pos_marker()
	}, bg="transparent")

	# click ion image
	observe({
		validate(need(input$ionimage_click, "Invalid click"))
		pos <- c(input$ionimage_click$x, input$ionimage_click$y)
		sv$xy(round(pos))
	})

	# click spectrum
	observe({
		validate(need(input$spectrum_click, "Invalid click"))
		sv$mz(input$spectrum_click$x)
	})

	# double-click ion image
	observe({
		validate(need(input$ionimage_dblclick, "Invalid double-click"))
		sv$ionimage_xylim(sv$xy_range())
	})

	# double-click spectrum
	observe({
		validate(need(input$spectrum_dblclick, "Invalid double-click"))
		sv$spectrum_massrange(sv$mz_range())
	})

	# brush ion image
	observe({
		validate(need(input$ionimage_brush, "Invalid brush"))
		p <- input$ionimage_brush
		xylim <- c(p$xmin, p$xmax, p$ymin, p$ymax)
		pos <- c((p$xmin + p$xmax) / 2, (p$ymin + p$ymax) / 2)
		sv$ionimage_xylim(xylim)
	})

	# brush spectrum
	observe({
		validate(need(input$spectrum_brush, "Invalid brush"))
		p <- input$spectrum_brush
		massrange <- c(p$xmin, p$xmax)
		mz <- (p$xmin + p$xmax) / 2
		sv$spectrum_massrange(massrange)
	})

	# change x/y image axes
	observe({
		xylim <- c(range(coord(data())[,sv$xy_names()[1]]),
			range(coord(data())[,sv$xy_names()[2]]))
		sv$ionimage_xylim(xylim)
	})

	return(sv)

}
