function loadCausalFlow(mutationName, regulatorName, biclusterName, phenotypeName) {
	// set the dimensions and margins of the graph
	const margin = {
		top: 40,
		right: 40,
		bottom: 40,
		left: 40
	}
	const boxWidth = 400 - margin.left - margin.right
	const boxHeight = 360 - margin.top - margin.bottom

	const width = 1170
	const height = 750

	const legendTextHeight = 20

	d3.select("#causal-flow svg").remove()

	document.getElementById("save-causal-flow").addEventListener("click", event => {
		saveSvgAsPng(document.getElementById("causal-flow-svg"), `${mutationName}-${regulatorName}-${biclusterName}.png`)
		saveSvg(document.getElementById("causal-flow-svg"), `${mutationName}-${regulatorName}-${biclusterName}.svg`)
	})

	// append the svg object to the body of the page
	const main = d3.select("#causal-flow")
		.append("svg")
		.attr("id", "causal-flow-svg")
		.attr("width", width)
		.attr("height", height)
		.style("background-color", "#FFF")

	let xAxisCount = 0, yAxisCount = 0
	const tooltip = d3.select("#causal-flow-tooltip")

	d3.json(`/bicluster-causal-analysis/${mutationName}/${regulatorName}/${biclusterName}/${phenotypeName}`, (json) => {
		main.attr("width", width)
			.attr("height", 750 + 40 + Math.ceil(json.scatter.legend.length / 2) * legendTextHeight) // adjust SVG height based off of legend
		
		// scatterplot
		const drawScatterplot = (left, top) => {
			const svg = main.append("g")
				.attr("transform", `translate(${margin.left + left}, ${margin.top + top})`)
			
			// stats text
			svg.append("text")
				.attr("x", boxWidth / 2)
				.attr("y", -5)
				.style("text-anchor", "middle")
				.style("font-family", "Arial")
				.style("font-size", "12pt")
				.text(`R = ${json.scatter.stats[0].toPrecision(3)}, p = ${json.scatter.stats[1].toExponential(2)}`)
			
			// x-axis text
			svg.append("text")
				.attr("x", boxWidth / 2)
				.attr("y", boxHeight + 45)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.text(`${regulatorName} Expression`)
			
			// y-axis text
			svg.append("text")
				.attr("x", 0)
				.attr("y", 0)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.style("transform", `rotate(270deg) translate(-${boxHeight / 2}px, -55px)`)
				.text(`${biclusterName} Expression`)
			
			// find min/max on each axis
			let minX = Number.MAX_VALUE, maxX = Number.MIN_VALUE, minY = Number.MAX_VALUE, maxY = Number.MIN_VALUE
			for (const datum of json.scatter.data) {
				for (const point of datum.data) {
					minX = Math.min(minX, point.x)
					maxX = Math.max(maxX, point.x)
					minY = Math.min(minY, point.y)
					maxY = Math.max(maxY, point.y)
				}
			}

			// x-axis
			const xMargin = (maxX - minX) * 0.05
			const x = d3.scaleLinear()
				.domain([minX - xMargin, maxX + xMargin])
				.range([0, boxWidth])

			svg.append("g")
				.attr("transform", "translate(0," + boxHeight + ")")
				.attr("class", `xAxis${xAxisCount}`)
				.call(d3.axisBottom(x).ticks(5))
			
			// x-axis styling
			d3.select(`.xAxis${xAxisCount} path`)
				.style("stroke-width", "0px")
			
			d3.selectAll(`.xAxis${xAxisCount} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.xAxis${xAxisCount} .tick text`)
				.style("font-size", "10pt")

			// y-axis
			const yMargin = (maxY - minY) * 0.05
			const y = d3.scaleLinear()
				.domain([minY - yMargin, maxY + yMargin])
				.range([boxHeight, 0])

			svg.append("g")
				.attr("class", `yAxis${yAxisCount}`)
				.call(d3.axisLeft(y).ticks(8))

			// y-axis styling
			d3.select(`.yAxis${yAxisCount} path`)
				.style("stroke-width", "0px")
			
			d3.selectAll(`.yAxis${yAxisCount} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.yAxis${yAxisCount} .tick text`)
				.style("font-size", "10pt")

			// add line
			svg.append("line")
				.style("stroke", "red")
				.style("stroke-width", 2)
				.attr("x1", x(json.scatter.regression[0][0]))
				.attr("y1", y(json.scatter.regression[0][1]))
				.attr("x2", x(json.scatter.regression[1][0]))
				.attr("y2", y(json.scatter.regression[1][1]))

			// add dots for mutation
			const nodes = svg.append("g")
				.selectAll("dot")
				.data(json.scatter.data[0].data)
				.enter()
				.append("circle")
				.attr("cx", d => x(d.x))
				.attr("cy", d => y(d.y))
				.attr("r", 4)
				.style("fill", d => d.color)
				.on("mouseover", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("r", 7)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0.9)
					
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
					
					if(d.z !== undefined) {
						tooltip.html(`<b>${d.name} (Mutated)</b><br />${regulatorName} Expression: ${d.x.toFixed(4)}<br />${biclusterName} Expression: ${d.y.toFixed(4)}<br />${d.phenotype}: ${d.z.toFixed(4)}`)
					}
					else {
						tooltip.html(`<b>${d.name} ${d.phenotype} (Mutated)</b><br />${regulatorName} Expression: ${d.x.toFixed(4)}<br />${biclusterName} Expression: ${d.y.toFixed(4)}`)
					}
					
					nodes.sort((a, b) => {
						return a === d ? 1 : -1
					})
				})
				.on("mouseout", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("r", 4)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0)
					
					tooltip.style("left", `-1000px`)
						.style("top", `-1000px`)
				})
			
			nodes.sort((a, b) => -1)
			
			// add crosses for wild type
			svg.append("g")
				.selectAll("dot")
				.data(json.scatter.data[1].data)
				.enter()
				.append("path")
				.attr("d", "M4 4 L-4 -4 M-4 4 L4 -4 Z")
				.attr("stroke", d => d.color)
				.attr("stroke-width", 3)
				.attr("transform", d => `translate(${x(d.x)}, ${y(d.y)})`)
				.on("mouseover", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("d", "M6 6 L-6 -6 M-6 6 L6 -6 Z")
						.attr("stroke-width", 4)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0.9)
					
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
					
					if(d.z !== undefined) {
						tooltip.html(`<b>${d.name} (Wild Type)</b><br />${regulatorName} Expression: ${d.x.toFixed(4)}<br />${biclusterName} Expression: ${d.y.toFixed(4)}<br />${d.phenotype}: ${d.z.toFixed(4)}`)
					}
					else {
						tooltip.html(`<b>${d.name} ${d.phenotype} (Wild Type)</b><br />${regulatorName} Expression: ${d.x.toFixed(4)}<br />${biclusterName} Expression: ${d.y.toFixed(4)}`)
					}
				})
				.on("mouseout", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("d", "M4 4 L-4 -4 M-4 4 L4 -4 Z")
						.attr("stroke-width", 3)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0)
					
					tooltip.style("left", `-1000px`)
						.style("top", `-1000px`)
				})
			
			// draw border around graph
			svg.append("rect")
				.attr("x", 0)
				.attr("y", 0)
				.attr("width", boxWidth + 1)
				.attr("height", boxHeight + 1)
				.style("stroke", "#000")
				.style("fill", "none")
				.style("stroke-width", "2px")

			// draw legend
			const colorDomain = d3.scaleLinear().range(
				json.scatter.legend.map(
					value => value.color
				)
			).domain(
				json.scatter.legend.map(
					value => json.scatter.legend.indexOf(value) + 1
				)
			)

			const legend = main.append("g")
				.attr("transform", `translate(${margin.left + left}, ${margin.top + top + boxHeight + 60})`)
				.attr("width", boxWidth)
				.attr("height", Math.ceil(json.scatter.legend.length / 2) * legendTextHeight)
				.selectAll("g")
				.data(colorDomain.domain().slice())
				.enter()
				.append("g")
				.attr("transform", (d, i) => `translate(${(i % 2) * boxWidth / 2}, ${Math.floor(i / 2) * legendTextHeight})`)
			
			legend.append("circle")
				.attr("r", 5)
				.style("fill", colorDomain)
				.attr("transform", (d, i) => "translate(10, 10)")

			legend.append("text")
				.data(json.scatter.legend.map(value => value.name))
				.attr("x", 20)
				.attr("y", 10)
				.attr("dy", ".35em")
				.text(d => d)
			
			xAxisCount++
			yAxisCount++
		}

		const drawBoxplot = (data, left, top, xAxisText, yAxisText) => {
			const svg = main.append("g")
				.attr("transform", `translate(${margin.left + left}, ${margin.top + top})`)

			// stats text
			svg.append("text")
				.attr("x", boxWidth / 2)
				.attr("y", -5)
				.style("text-anchor", "middle")
				.style("font-family", "Arial")
				.style("font-size", "12pt")
				.text(`T = ${data.stats[0].toPrecision(3)}, p = ${data.stats[1].toExponential(2)}`)

			// x-axis text
			svg.append("text")
				.attr("x", boxWidth / 2)
				.attr("y", boxHeight + 45)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.text(xAxisText)

			if(yAxisText instanceof Array) {
				// y-axis text
				for(let i = 0; i < yAxisText.length; i++) {
					const text = yAxisText[i]
					svg.append("text")
						.attr("x", 0)
						.attr("y", 0)
						.style("text-anchor", "middle")
						.style("font-weight", "bold")
						.style("font-family", "Arial")
						.style("font-size", "14pt")
						.style("transform", `rotate(270deg) translate(-${boxHeight / 2}px, -${55 + (yAxisText.length - i - 1) * 25}px)`)
						.text(text)
				}
			}
			else {
				// y-axis text
				svg.append("text")
					.attr("x", 0)
					.attr("y", 0)
					.style("text-anchor", "middle")
					.style("font-weight", "bold")
					.style("font-family", "Arial")
					.style("font-size", "14pt")
					.style("transform", `rotate(270deg) translate(-${boxHeight / 2}px, -55px)`)
					.text(yAxisText)
			}
			
			// x-axis
			const x = d3.scaleBand()
				.range([0, boxWidth])
				.domain(["WT", "Mutated"])
				.paddingInner(1)
				.paddingOuter(0.5)
			
			svg.append("g")
				.attr("class", `xAxis${xAxisCount}`)
				.attr("transform", "translate(0," + boxHeight + ")")
				.call(d3.axisBottom(x))
			
			// x-axis styling
			d3.select(`.xAxis${xAxisCount} path`)
				.style("stroke", "transparent")
			
			d3.selectAll(`.xAxis${xAxisCount} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.xAxis${xAxisCount} .tick text`)
				.style("font-size", "11pt")
			
			let minY = Number.MAX_VALUE, maxY = Number.MIN_VALUE
			for(const datum of data.data) {
				minY = Math.min(minY, datum.low)
				maxY = Math.max(maxY, datum.high)
			}

			// y-axis
			const yMargin = (maxY - minY) * 0.05
			const y = d3.scaleLinear()
				.domain([minY - yMargin, maxY + yMargin])
				.range([boxHeight, 0])

			svg.append("g")
				.attr("class", `yAxis${yAxisCount}`)
				.call(d3.axisLeft(y).ticks(8))

			// y-axis styling
			d3.select(`.yAxis${yAxisCount} path`)
				.style("stroke-width", "0px")
			
			d3.selectAll(`.yAxis${yAxisCount} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.yAxis${yAxisCount} .tick text`)
				.style("font-size", "10pt")

			// enter data
			svg.selectAll("boxplot")
				.data(data.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name)))
				.attr("x2", d => Math.floor(x(d.name)))
				.attr("y1", d => Math.floor(y(d.low)))
				.attr("y2", d => Math.floor(y(d.high)))
				.attr("stroke", "black")
				.attr("stroke-width", "2px")

			const length = 50
			// max horizontal line
			svg.selectAll("boxplot")
				.data(data.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.high)))
				.attr("y2", d => Math.floor(y(d.high)))
				.attr("stroke", "black")
				.attr("stroke-width", "2px")
			
			// min horizontal line
			svg.selectAll("boxplot")
				.data(data.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.low)))
				.attr("y2", d => Math.floor(y(d.low)))
				.attr("stroke", "black")
				.attr("stroke-width", "2px")
			
			// drawing the boxes
			svg.selectAll("boxplot")
				.data(data.data)
				.enter()
				.append("rect")
				.attr("x", d => Math.floor(x(d.name) - length))
				.attr("y", d => Math.floor(y(d.q3)))
				.attr("height", d => Math.floor(y(d.q1) - y(d.q3)))
				.attr("width", length * 2)
				.attr("stroke", "black")
				.attr("stroke-width", "2px")
				.style("fill", d => d.fillColor)
				.on("mouseover", d => {
					tooltip.transition()
						.duration(200)
						.style("opacity", 0.9)
					
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
					
					tooltip.html(`<b>${d.name}</b><br />Maximum: ${d.high.toFixed(4)}<br />Upper Quartile: ${d.q3.toFixed(4)}<br />Median: ${d.median.toFixed(4)}<br />Loewr Quartile: ${d.q1.toFixed(4)}<br />Minimum: ${d.low.toFixed(4)}`)
				})
				.on("mousemove", d => {
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
				})
				.on("mouseout", d => {
					tooltip.transition()
						.duration(200)
						.style("opacity", 0)
					
					tooltip.style("left", `-1000px`)
						.style("top", `-1000px`)
				})
			
			// median line
			svg.selectAll("boxplot")
				.data(data.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.median)))
				.attr("y2", d => Math.floor(y(d.median)))
				.attr("stroke", "black")
				.attr("stroke-width", "3px")
			
			// draw border around graph
			svg.append("rect")
				.attr("x", 0)
				.attr("y", 0)
				.attr("width", boxWidth + 1)
				.attr("height", boxHeight + 1)
				.style("stroke", "#000")
				.style("fill", "none")
				.style("stroke-width", "2px")
			
			xAxisCount++
			yAxisCount++
		}

		drawBoxplot(json.mutation_gene_expression, 40, 0, `${regulatorName} Mutational Status`, `${regulatorName} Expression`)
		drawBoxplot(json.bicluster_eigengene_expression, 800, 0, `${regulatorName} Mutational Status`, `${biclusterName} Expression`)
		drawBoxplot(json.residual, 800, 400, `${regulatorName} Mutational Status`, [`${biclusterName} Expression`, `Conditioned On ${regulatorName} Expression`])
		drawScatterplot(40, 400)

		// handle the middle decorations
		{
			const svg = main.append("g")
			const top = 260
			
			// mutation
			svg.append("path")
				.attr("d", "M0 0 L40 -90 L0 -65 L-40 -90 L0 0 Z")
				.attr("fill", "#C92CD4")
				.attr("stroke", "#000")
				.attr("stroke-width", "2px")
				.style("transform", `translate(${width / 2}px, ${top}px)`)
			
			svg.append("text")
				.attr("x", width / 2)
				.attr("y", 30 + top)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.text(mutationName)

			// regulator
			svg.append("path")
				.attr("d", "M0 0 L40 80 L-40 80 L0 0 Z")
				.attr("fill", "#FFAA00")
				.attr("stroke", "#000")
				.attr("stroke-width", "2px")
				.style("transform", `translate(${width / 2}px, ${90 + top}px)`)
			
			svg.append("text")
				.attr("x", width / 2)
				.attr("y", 195 + top)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.text(regulatorName)
			
			// bicluster
			const biclusterWidth = 80
			svg.append("rect")
				.attr("x", width / 2 - biclusterWidth / 2)
				.attr("y", 254 + top)
				.attr("width", biclusterWidth)
				.attr("height", biclusterWidth)
				.style("stroke", "#000")
				.style("fill", "#CCC")
				.style("stroke-width", "2px")
			
			svg.append("text")
				.attr("x", width / 2)
				.attr("y", 359 + top)
				.style("text-anchor", "middle")
				.style("font-weight", "bold")
				.style("font-family", "Arial")
				.style("font-size", "14pt")
				.text(biclusterName)
			
			// right line + arrow
			svg.append("line")
				.attr("x1", width / 2 + 100)
				.attr("y1", top - 90)
				.attr("x2", width / 2 + 85)
				.attr("y2", top - 90)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 + 100)
				.attr("y1", top + 360)
				.attr("x2", width / 2 + 85)
				.attr("y2", top + 360)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 + 100)
				.attr("y1", top - 90 - 2)
				.attr("x2", width / 2 + 100)
				.attr("y2", top + 360 + 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 + 100)
				.attr("y1", top + 135)
				.attr("x2", width / 2 + 125)
				.attr("y2", top + 135)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 + 125)
				.attr("y1", top - 75 - 2)
				.attr("x2", width / 2 + 125)
				.attr("y2", top + 135 + 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 + 125)
				.attr("y1", top - 75)
				.attr("x2", width / 2 + 150)
				.attr("y2", top - 75)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("path")
				.attr("d", "M0 0 L-18 -9 L-18 9 L0 0 Z")
				.attr("fill", "#000")
				.style("transform", `translate(${width / 2 + 168}px, ${top - 75}px)`)
			
			// arrow from mutation to regulator
			svg.append("line")
				.attr("x1", width / 2)
				.attr("y1", top + 40)
				.attr("x2", width / 2)
				.attr("y2", top + 71)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("path")
				.attr("d", "M0 0 L-9 -18 L9 -18 L0 0 Z")
				.attr("fill", "#000")
				.style("transform", `translate(${width / 2}px, ${top + 89}px)`)
			
			
			// arrow from regulator to bicluster
			svg.append("line")
				.attr("x1", width / 2)
				.attr("y1", top + 205)
				.attr("x2", width / 2)
				.attr("y2", top + 236)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")

			svg.append("path")
				.attr("d", "M0 0 L-9 -18 L9 -18 L0 0 Z")
				.attr("fill", "#000")
				.style("transform", `translate(${width / 2}px, ${top + 254}px)`)

			// left line
			svg.append("line")
				.attr("x1", width / 2 - 85)
				.attr("y1", top - 90)
				.attr("x2", width / 2 - 100)
				.attr("y2", top - 90)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 85)
				.attr("y1", top + 195)
				.attr("x2", width / 2 - 100)
				.attr("y2", top + 195)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100)
				.attr("y1", top - 90 - 2)
				.attr("x2", width / 2 - 100)
				.attr("y2", top + 195 + 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100)
				.attr("y1", top + 50)
				.attr("x2", width / 2 - 125)
				.attr("y2", top + 50)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 125)
				.attr("y1", top + 50 + 2)
				.attr("x2", width / 2 - 125)
				.attr("y2", top - 75 - 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 125)
				.attr("y1", top - 75)
				.attr("x2", width / 2 - 150)
				.attr("y2", top - 75)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")

			svg.append("path")
				.attr("d", "M0 0 L18 9 L18 -9 L0 0 Z")
				.attr("fill", "#000")
				.style("transform", `translate(${width / 2 - 168}px, ${top - 75}px)`)
			
			// another line that goes somewhere
			const lineLength = 17
			svg.append("line")
				.attr("x1", width / 2 - 100)
				.attr("y1", top + 90)
				.attr("x2", width / 2 - 100 - lineLength)
				.attr("y2", top + 90)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100)
				.attr("y1", top + 360)
				.attr("x2", width / 2 - 100 - lineLength)
				.attr("y2", top + 360)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100 - lineLength)
				.attr("y1", top + 90 - 2)
				.attr("x2", width / 2 - 100 - lineLength)
				.attr("y2", top + 360 + 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100 - lineLength)
				.attr("y1", top + 225)
				.attr("x2", width / 2 - 100 - lineLength * 2)
				.attr("y2", top + 225)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100 - lineLength * 2)
				.attr("y1", top + 225 - 2)
				.attr("x2", width / 2 - 100 - lineLength * 2)
				.attr("y2", top + 325 + 2)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("line")
				.attr("x1", width / 2 - 100 - lineLength * 2)
				.attr("y1", top + 325)
				.attr("x2", width / 2 - 100 - lineLength * 3)
				.attr("y2", top + 325)
				.attr("stroke", "black")
				.attr("stroke-width", "4px")
			
			svg.append("path")
				.attr("d", "M0 0 L18 9 L18 -9 L0 0 Z")
				.attr("fill", "#000")
				.style("transform", `translate(${width / 2 - 100 - lineLength * 3 - 18}px, ${top + 325}px)`)
		}
	})
}