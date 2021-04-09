function plotAchillesBoxplot(bicluster, containerWidth) {
	d3.select("#achilles-boxplot svg").remove()
	
	const main = d3.select("#achilles-boxplot")
		.append("svg")
		.attr("id", "achilles-boxplot-svg")
		.style("background-color", "#FFF")
	
	const boxWidth = containerWidth
	const boxHeight = 600
	const bottomPad = 160

	const legendTextHeight = 20

	const tooltip = d3.select("#achilles-tooltip")

	document.getElementById("save-achilles-boxplot").addEventListener("click", event => {
		saveSvgAsPng(document.getElementById("achilles-boxplot-svg"), `${bicluster}.png`)
		saveSvg(document.getElementById("achilles-boxplot-svg"), `${bicluster}.svg`)
	})
	
	d3.json(`/achilles-bicluster/${bicluster}`, (json) => {
		json.data = json.data.sort((a, b) => a.name.localeCompare(b.name))
		main.attr("width", boxWidth + 60)
			.attr("height", 10 + boxHeight + bottomPad)
		
		const svg = main.append("g")
			.attr("transform", `translate(50, 10)`)
		
		// x-axis
		const x = d3.scaleBand()
			.range([0, boxWidth])
			.domain(json.data.map(a => a.name))
			.paddingInner(1)
			.paddingOuter(0.5)
	
		svg.append("g")
			.attr("class", `xAxis${bicluster}`)
			.attr("transform", "translate(0," + boxHeight + ")")
			.call(d3.axisBottom(x))
		
		// x-axis styling
		d3.select(`.xAxis${bicluster} path`)
			.style("stroke", "transparent")
		
		d3.selectAll(`.xAxis${bicluster} .tick`)
			.style("stroke-width", "2px")
		
		d3.selectAll(`.xAxis${bicluster} .tick text`)
			.style("font-size", "11pt")
			.style("transform", "rotate(270deg) translate(-10px, -14px)")
			.style("text-anchor", "end")
		
		// y-axis
		const minimumScore = json.data.reduce((minimum, current) => current.low < minimum ? current.low : minimum, 500)
		const minY = Math.min(-2, minimumScore), maxY = 2
		const y = d3.scaleLinear()
			.domain([minY, maxY])
			.range([boxHeight, 0])
			.nice()

		svg.append("g")
			.attr("class", `yAxis${bicluster}`)
			.call(d3.axisLeft(y).ticks(8).tickSize(-boxWidth))

		// y-axis styling
		d3.select(`.yAxis${bicluster} path`)
			.style("stroke-width", "0px")
		
		d3.selectAll(`.yAxis${bicluster} line`)
			.style("stroke", "#AAA")
		
		d3.selectAll(`.yAxis${bicluster} .tick`)
			.attr("transform", d => `translate(0, ${Math.floor(y(d))})`)
			.style("stroke-width", "2px")
		
		d3.selectAll(`.yAxis${bicluster} .tick text`)
			.style("font-size", "10pt")
		
		const tooltipify = (o) => {
			return o.on("mouseover", d => {
				tooltip.transition()
					.duration(200)
					.style("opacity", 0.9)
				
				tooltip.style("left", `${d3.event.clientX}px`)
					.style("top", `${d3.event.clientY}px`)
				
				tooltip.html(`<b>${d.name}</b><br />Maximum: ${d.high.toFixed(4)}<br />Upper Quartile: ${d.q3.toFixed(4)}<br />Median: ${d.median.toFixed(4)}<br />Lower Quartile: ${d.q1.toFixed(4)}<br />Minimum: ${d.low.toFixed(4)}`)
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
		}

		const length = Math.max((containerWidth / json.data.length - 20) / 2, 5)
		// line through -1
		svg.append("line")
			.attr("x1", 0)
			.attr("x2", boxWidth)
			.attr("y1", Math.floor(y(-1)))
			.attr("y2", Math.floor(y(-1)))
			.attr("stroke", "red")
			.attr("stroke-width", "2px")
		
		// line through 0
		svg.append("line")
			.attr("x1", 0)
			.attr("x2", boxWidth)
			.attr("y1", Math.floor(y(0)))
			.attr("y2", Math.floor(y(0)))
			.attr("stroke", "#555")
			.attr("stroke-width", "2px")
		
		// y-axis text
		svg.append("text")
			.attr("x", 0)
			.attr("y", 0)
			.style("text-anchor", "middle")
			.style("font-weight", "bold")
			.style("font-family", "Arial")
			.style("font-size", "14pt")
			.style("transform", `rotate(270deg) translate(${-y(0)}px, -35px)`)
			.text(`Gene Effect Ceres`)

		const boxplotOutlineColor = "#333333"

		// enter data
		tooltipify(
			svg.selectAll("boxplot")
				.data(json.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name)))
				.attr("x2", d => Math.floor(x(d.name)))
				.attr("y1", d => Math.floor(y(d.low)))
				.attr("y2", d => Math.floor(y(d.high)))
				.attr("stroke", boxplotOutlineColor)
				.attr("stroke-width", "2px")
		)

		// max horizontal line
		tooltipify(
			svg.selectAll("boxplot")
				.data(json.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.high)))
				.attr("y2", d => Math.floor(y(d.high)))
				.attr("stroke", boxplotOutlineColor)
				.attr("stroke-width", "2px")
		)
		
		// min horizontal line
		tooltipify(
			svg.selectAll("boxplot")
				.data(json.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.low)))
				.attr("y2", d => Math.floor(y(d.low)))
				.attr("stroke", boxplotOutlineColor)
				.attr("stroke-width", "2px")
		)
		
		// drawing the boxes
		tooltipify(
			svg.selectAll("boxplot")
				.data(json.data)
				.enter()
				.append("rect")
				.attr("x", d => Math.floor(x(d.name) - length))
				.attr("y", d => Math.floor(y(d.q3)))
				.attr("height", d => Math.floor(y(d.q1) - y(d.q3)))
				.attr("width", length * 2)
				.attr("stroke", boxplotOutlineColor)
				.attr("stroke-width", "2px")
				.style("fill", "#FFFFFF")
		)
		
		// median line
		tooltipify(
			svg.selectAll("boxplot")
				.data(json.data)
				.enter()
				.append("line")
				.attr("x1", d => Math.floor(x(d.name) - length))
				.attr("x2", d => Math.floor(x(d.name) + length))
				.attr("y1", d => Math.floor(y(d.median)))
				.attr("y2", d => Math.floor(y(d.median)))
				.attr("stroke", boxplotOutlineColor)
				.attr("stroke-width", "2px")
		)

		// dots for each cell line
		for(const datum of json.data) {
			const gene = datum.name
			const nodes = svg.selectAll("dot")
				.data(datum.data)
				.enter()
				.append("circle")
				.attr("cx", d => Math.floor(x(gene)))
				.attr("cy", d => Math.floor(y(d.score)))
				.attr("r", 3)
				.style("fill", d => d.color)
				.on("mouseover", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("r", 4.5)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0.9)
					
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
					
					tooltip.html(`<b>${d.cell_line} (${d.subtype})</b><br />Effect: ${d.score.toFixed(4)}<br />Dependency: ${d.p_value.toFixed(4)}`)
					
					nodes.sort((a, b) => {
						return a === d ? 1 : -1
					})
				})
				.on("mouseout", function(d) {
					d3.select(this)
						.transition()
						.duration(200)
						.attr("r", 3)
					
					tooltip.transition()
						.duration(200)
						.style("opacity", 0)
					
					tooltip.style("left", `-1000px`)
						.style("top", `-1000px`)
				})
			
			nodes.sort((a, b) => -1)
		}

		// draw legend
		const colorDomain = d3.scaleLinear().range(
			json.legend.map(
				value => value.color
			)
		).domain(
			json.legend.map(
				value => json.legend.indexOf(value) + 1
			)
		)

		const legend = main.append("g")
			.attr("transform", `translate(50, ${boxHeight + 120})`)
			.attr("width", 300)
			.attr("height", Math.ceil(json.legend.length / 2) * legendTextHeight)
			.selectAll("g")
			.data(colorDomain.domain().slice())
			.enter()
			.append("g")
			.attr("transform", (d, i) => `translate(${(i % 2) * 300 / 2}, ${Math.floor(i / 2) * legendTextHeight})`)
		
		legend.append("circle")
			.attr("r", 5)
			.style("fill", colorDomain)
			.attr("transform", (d, i) => "translate(10, 10)")
		
		legend.append("text")
			.data(json.legend.map(value => value.subtype))
			.attr("x", 20)
			.attr("y", 10)
			.attr("dy", ".35em")
			.text(d => d)
	})
}