<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">

  <title>Interactive heatmap</title>
  <meta name="description" content="Interactive heatmap and dot plot">
  <meta name="author" content="Dorin-Mirel Popescu and Issac Goh">
  <style type = 'text/css'>
  	::-webkit-scrollbar {
    	-webkit-appearance: none;
    	width: 2px;
	}
	::-webkit-scrollbar-thumb {
    	border-radius: 1px;
    	background-color: rgba(0,0,0,.5);
    	box-shadow: 0 0 1px rgba(255,255,255,.5);
	}
  	#cell-menu {
    	overflow-x: auto;
    	white-space: nowrap;
    	margin: .5em;
    	padding: .5em;
  	}
  	.menuContainer {
  		border-width: thin;
  		border-style: solid;
  		margin: 4px;
  		padding: 10px;
  	}
  	#instructionBtn {
  		margin-bottom: 1em;
  		font-size: 1em;
  		background-color: #00bfff;
    	border: none;
    	color: white;
    	padding: 10px 10px;
    	text-align: center;
    	text-decoration: none;
    	cursor: pointer;
  	}
  	#end_div {
  		margin-bottom: 3em;
  	}
  	#instructions_div {
  		margin:0; padding:0;
  	}
  	body {
  font-family: Avenir, Arial, sans-serif;
}
  </style>
</head>

<body>
	<input id = 'instructionBtn' type = 'button' value = 'Show description and instructions' onclick = 'showHideInstructions(value)' />
	<div id = 'instructionBoard'>
		<div id = 'dataset_div'><b>Data set: </b></div>
		<div><b>Description: </b>An interactive environment for exploring the expression of multiple genes in a data set using dot plots and heatmaps.</div>
		<div>
			<b>Instructions: </b>
			<div id = 'instructions_div'><ul>
				<li>Switch between cell types on x axis and gene names on x axis using the drop-down menu at "Layout"</li>
				<li>Switch between dot plot and heatmap using the drop down menu at "Plot type"</li>
				<li>The "Cell type selection menu" allows switching on/off of particular cell types. Changes will be reflected in the plots immediately. If there are too many cell types, use to horizontal scrollbar to navigate the list of cell types.</li>
				<li>In the "Cell type selection menu" there is also the possibility to write or copy/paste the list desired cell types</li>
				<li>The "Gene selection menu" contains functionalities for editing the list of genes that are plotted</li>
				<li>Each gene being plotted has a dedicated field which indicates its name. As your mouse moves over a gene its name is boldened and a "remove" button will appear, allowing you yo remove genes from the plot. You can also clear the list using the "Remove All" button.</li>
				<li>Using the dropdown list you can choose to search for Gene symbols, search for diseases, or type/paste a list of gene symbols.</li>
				<li>As you search for gene symbols they will be added to the list. To add a gene, click its name in the search results. Alternatively, if you've typed a full gene symbol pressing enter will also add it to the list.</li>
				<li>If you search for diseases, the number of related genes is shown in brackets after the gene name. Clicking an option in the list will remove any currently selected genes and add all the disease-related genes to the plot.</li>
				<li>Finally, a list of genes can be pasted or typed in the text box. The list should have the gene names separated by commas.</li>
				<li><b>You can click-chose any gene or cell type directly on the plot and reorder them by dragging</b></li>
				<li>Plots can be saved as png files by right clicking on plot area and choosing "Save as"</li>
				<li>At the bottom of the page there is a list of genes and a list of cell types used in the plot. These lists can be copy/pasted and saved for easy restoring of the plots in a future session.</li>
			</ul></div>
		</div>
		<div><b>Applications: </b>FACS sorting panel design; data annotation; cell type validation; exploring groups of genes relevant to a function (e.g. cycling genes)</div>
		<div id = 'end_div'><hr/></div>
	</div>

    <form>
	<div>Layout: <select onchange = 'setLayout(this.value)'>
		<option value = 'xcells'>'Cells on X; Genes on Y'<option>
		<option value = 'xgenes' selected>'Cells on Y; Genes on X'<option>
	</select>
	Plot type: <select onchange = 'setPlotType(this.value)'>
    <option value = 'spotplot'>Dot plot</option>
    <option value = 'heatplot'>Heatmap</option>
	</select>
    
    
    Categorical variable <select id = "cat_var" name = "variable" onchange= "this.form.submit()">
    <option value="" selected>Selected variable--></option>
      <?php
          $data_files = scandir('./main');
          foreach($data_files as $data_file){
              $data_file_name = str_replace(".html","",$data_file);
              $data_file_name = str_replace(".php","",$data_file_name);
              $data_file_name = str_replace(".csv","",$data_file_name);
              $data_file_name = str_replace(".txt","",$data_file_name);
              if (($data_file_name !== '.') && ($data_file_name !=='..')) {
              $filename = basename($data_file_name);
              echo "<option value='". $filename ."'>". $filename ."</option>";
          }
        }
    ?>
<?php
    if(isset($_GET["variable"])){
        $variable=$_GET["variable"];
        $variable_name = $variable . '.txt';
        $data_var = file_get_contents('./main/'.$variable_name, FILE_USE_INCLUDE_PATH);
    }else{
        $variable = scandir('./main')[2];
        $data_var = file_get_contents('./main/'.$variable, FILE_USE_INCLUDE_PATH);
        
    }
    ?>
    </select>
<label for="cat_var">Variable selected: "<?php echo $variable?>"
</label>
    </form>
    </div>



        
	<div class = 'menuContainer'>
		<b>Cell type selection menu</b><hr>
		&nbsp;&nbsp;Toggle on/off each cell type individually:
		<div id ='cell-menu'></div>
		&nbsp;&nbsp;Or enter here a list of cell types here and press enter: 
		<input type ='text' id = 'cellTypeListReader' size = '20' />
	</div>
	<div class = 'menuContainer'><div><b>Gene selection menu</b></div><hr>
	    <gene-selector multiple autoselect buttons="remove" data-init="CD34, SPINK2, PRSS57, ALAS2, GYPA, KLF1, GP9, PLEK, ITGA2B, VPREB1, IGLL1, CD79A, CLEC10A, CD1C, HLA-DPA1, GATA2, HDC, PRG2, NCF1, ITGAM, PGLYRP1, CD14, CD68, CD52, CD3E, TRBC1, GZMA, MMP9, KDR, CTHRC1">
		<endpoint name="genes" function="findGenes" mode="list" property="results" empty="msg:No matching genes for this dataset.">Search Gene Symbols</endpoint>
		<endpoint name="diseases" url="/datasets/fbm_updated/diseases/fetch_disease_data.php" preprocessor="verifyDiseaseData" parameter="term" property="data" mode="groups" action="replace" empty="msg:No matching diseases in this database." searchblank>Search for Diseases</endpoint>
		<endpoint name="list" function="processGeneList" mode="list" property="results" empty="noaction" trigger="enter">Paste a list of genes</endpoint>
	    </gene-selector>
	</div>
	<div><canvas id = 'interactive-canvas' width = '700' height = '300' ></canvas></div>
	<div>
		<div>List of cell types: 
			<ul><li id = 'output_div_cells' ></li></ul>
		</div>
		<div>List of genes: 
			<ul><li id = 'output_div_genes' ></li></ul>
		</div>
	</div>
<script type="module">
import { ValueSelector } from '/datasets/fbm_updated/diseases/ValueSelector.js';
customElements.define("gene-selector", ValueSelector);
</script>

<script>
function findGenes(term) {
	return {
		results: Object.keys(expression_data).filter(k => k.startsWith(term))
	}
}

function processGeneList(term) {
    let genes = term.split(",");
    genes = genes.map(s => s.trim());
    let selector = document.querySelector("gene-selector");
    selector.find(".header").textContent = "Your saved gene list:";
    selector.find(".holder").classList.add("hidden");
    selector.removeAll();
    selector.push(genes);
    selector.find("input[type='text']").value = "";
    return { results: [] } // returning an empty list for this endpoint will prevent any further action.
}

function verifyDiseaseData(data) {
	let diseases = Object.keys(data.data);
	diseases.forEach(d => { 
	    data.data[d] = verifyData(data.data[d]);
	    if (data.data[d].length == 0) delete data.data[d];
	});
	return data;
}

function verifyData(data) {
	return data.filter(g => Boolean(g) && Boolean(expression_data[g]));
}
</script>
	<script type = 'text/javascript'>
		// global parameters
		var canvas             = document.getElementById('interactive-canvas'),	
			context            = canvas.getContext('2d'),
			layout             = 'xgenes', //other value is xcells
			plotType           = 'spotplot', // other value is heatmap
			selectedCellTypes  = [], // a named array. names are cell names, and values are index position on the axis
			selectedGenes      = [], // a named array. names are gene names, and values are index position on the axis
			row_coordinates    = [], // placeholder for row positions
			col_coordinates    = [], // placeholder for column positions
			padX               = null, // number of pixels at the left of the plot to fit the longest ylabel
			padY               = null, // number of pixels at bottom of the plot to fit the longest xlabel
			canvasW            = 700,
			canvasH            = 400,
			cellMenu           = document.getElementById('cell-menu'),
			geneMenu           = document.getElementById('gene-menu'),
			geneListReader     = document.getElementById('geneListReader'),
			cellTypeListReader = document.getElementById('cellTypeListReader'),
			instructionBtn     = document.getElementById('instructionBtn'),
			instructionBoard   = document.getElementById('instructionBoard'),
			dataset_div        = document.getElementById('dataset_div'),
			output_div_cells   = document.getElementById('output_div_cells'),
			output_div_genes   = document.getElementById('output_div_genes'),
			partyStarted       = false, // set to true when document loaded and the procedure have been put in placed
			tileDim            = 30,
			labelsFont         = '20px Avenir', // font for wrinting plot labels
			colorStart         = '#98f5ff', // color-start for the heatmap color gradient
			colorMiddle        = '#D2691E', // color-middle for the heatmap color gradient
			colorStop          = '#8b008b', // color-stop for the heatmap color gradoent
			mouseXStart        = null, // placeholder for the X position where a mouse click triggered an event
			mouseYStart        = null, // placeholder for the Y position where a mouse click triggered an event
			mouseXCurrent      = null, // placeholder for the cursor X position during  a drag event
			mouseYCurrent      = null, // placeholder for the cursor Y position during  a drag event
			choosenRowIndex    = null, // placeholder to store index of row during dragging
			choosenColIndex    = null, // placeholder to store index of column during dragging
			selectedRowPos     = null, // placeholder for the position of the selected row. Used for updating its corresponding real row position during dragging
			selectedColPos     = null; // placeholder for the position of the selected col. Used for updating its corresponding real col position during dragging
		
		// function to convert hex color string to rgb values
		function hexColorStringToRGB(hex_string){
			var r = parseInt('0x' + hex_string.slice(1, 3)),
				g = parseInt('0x' + hex_string.slice(3, 5)),
				b = parseInt('0x' + hex_string.slice(5, 7));
			return {'r':r, 'g':g, "b":b}
		}
		
		// conver the heatmap colors from hex strings to rgb values. Values are to be usedfor interpolation when computing the colors for heatmap
		colorStart  = hexColorStringToRGB(colorStart);
		colorMiddle = hexColorStringToRGB(colorMiddle);
		colorStop   = hexColorStringToRGB(colorStop);
			
		function lengthOfNamedArray(namedArray){
			var nal = 0;
			for (key in namedArray){nal++}
			return nal;
		}
			
		// data placeholders
		expression_data  = [] // an array of arrays (2D matrix) -> for each gene name insert one array with values for all the cell names
		cell_names_order = [] // array containing cell names in the same order as they appear in the above matrix
        
        var name = <?php echo($data_var); ?>;

		// put the name of the dataset
		dataset_div.innerHTML += dataset_name;
		
		// hide instructions/description
		instructionBoard.style.display = 'none'
		
		// compute the number of genes - useful later in the script
		// and get the gene nanes - use during the script for data handling
		var noGenes = 0,
			gene_names = [];
		for (key in expression_data){noGenes += 1; gene_names.push(key)}
		
		// populate selectedCellTypes
		cell_names.forEach(function(cell_name, index){selectedCellTypes[cell_name] = index})
		
		// function to show/hide description and instructions
		function showHideInstructions(val){
			if (val == 'Show description and instructions'){
				instructionBtn.value = 'Hide description and instructions'
				instructionBoard.style.display = 'block'
			}else{
				instructionBtn.value = 'Show description and instructions'
				instructionBoard.style.display = 'none'
			}
		}
		
		// function to allocate position for each row and each column
		// this function also resize the canvas to fit the number of columns and rows
		// when ever there is a re-allocation of columns and row the canvas should be redraw - therefore this function ends with a draw() call
		function allocateRowAndColumnsPositions(){
			// if layout is xcells, the the number of columns is the number of selected cells and the number of rows is the number of selected genes
			var rowLabels = [],
				colLabels = [];
			if (layout == 'xcells'){
				var ncols = lengthOfNamedArray(selectedCellTypes),
					nrows = lengthOfNamedArray(selectedGenes);
				for (key in selectedCellTypes){colLabels.push(key)}
				for (key in selectedGenes){rowLabels.push(key)}
			}else { // else, just the other way around
				var nrows = lengthOfNamedArray(selectedCellTypes),
					ncols = lengthOfNamedArray(selectedGenes);
				for (key in selectedCellTypes){rowLabels.push(key)}
				for (key in selectedGenes){colLabels.push(key)}
			}
			
			// if we have no cells or no genes clear the canvas and return.
			if (rowLabels.length == 0 || colLabels.length == 0) {
			    context.fillStyle = '#ffffff' // should change to white at the end
    			context.fillRect(0, 0, canvas.width, canvas.height)
    			return;
    	    }

			// update padX based on row labels
			context.font = labelsFont;
			padX = 5 + rowLabels.map(function(name){return context.measureText(name).width}).reduce(function(num1, num2){return Math.max(num1, num2)})
			// update padY based on column labels
			padY = 5 + 0.69 * colLabels.map(function(name){return context.measureText(name).width}).reduce(function(num1, num2){return Math.max(num1, num2)})
			// update canvas width and height
			canvasW       = padX + tileDim * ncols
			canvasH       = padY + tileDim * nrows
			canvas.width  = canvasW;
			canvas.height = canvasH
			context       = canvas.getContext('2d')
			// update the col and row coordinates
			col_coordinates = []
			row_coordinates = []
			for(i = 0; i < ncols; i++){col_coordinates.push(padX + tileDim * i)}
			for(i = 0; i < nrows; i++){row_coordinates.push(tileDim * i)}
			// and finally draw the canvas
			draw()
		}
			
		function setLayout(value){
			layout = value
			allocateRowAndColumnsPositions()
		}
		
		function setPlotType(value){
			plotType = value;
			draw()
		}
		
		function updateSelectedCellTypes(){
			// must handle addition and removal of cell types differently
			var newSelectedCellTypes = [],
				cell_buttons = cellMenu.children;
			for(i = 0; i < cell_buttons.length; i++){
				var spanID = cell_buttons[i];
				if (spanID.children[0].checked){
					newSelectedCellTypes[spanID.children[0].value] = 0
				}
			}
			if (lengthOfNamedArray(selectedCellTypes) > lengthOfNamedArray(newSelectedCellTypes)){ // condition for removal of cell type
				// must get the position of the removed cell name
				// then decrement all position above it by 1
				var disjointPosition = 0;
				for (key in selectedCellTypes){
					if (!(key in newSelectedCellTypes)){
						disjointPosition = selectedCellTypes[key]
					}
				}
				//now update the newCellSelectedCellTypes position
				for(key in newSelectedCellTypes){
					var position = selectedCellTypes[key]
					if (position > disjointPosition){
						newSelectedCellTypes[key] = position - 1
					}else{
						newSelectedCellTypes[key] = position
					}
				}
				selectedCellTypes = newSelectedCellTypes;
			}else if(lengthOfNamedArray(selectedCellTypes) < lengthOfNamedArray(newSelectedCellTypes)){
				// must get the additional cell name
				// then insert the additional cell name at the end of selectedCellTypes and assign it increment of the the highest Position
				var insertedCellType = null;
				for (key in newSelectedCellTypes){
					if (!(key in selectedCellTypes)){insertedCellType = key}
				}
				var highestPos = 0;
				for(key in selectedCellTypes){highestPos = Math.max(highestPos, selectedCellTypes[key])}
				selectedCellTypes[insertedCellType] = highestPos + 1;
			}
			// and the block of instruction sends with a call to allocateRowAndColumnsPositions:
			allocateRowAndColumnsPositions()
		}
		
		// this function inserts a random gene on the table. 
		// because it changes the workable data, at the end it calls the allocateRowAndColumnsPositions()
		// which in turn call draw()
		function generateNewGeneButton(){
			expressionVariance = 0
			while (expressionVariance < .3 * cell_names.length){
				var index          = parseInt(Math.random() * noGenes),
					gene_name      = gene_names[index]
					expressionMean = expression_data[gene_name].reduce(function(a, b){return a + b}) / cell_names.length
					diffs          = expression_data[gene_name].map(function(a){return Math.abs(a - expressionMean)}),
				expressionVariance = diffs.reduce(function(a, b){return a + b })
				// if gene already present in list, loop some more - otherwise there will be problems
				if (gene_name in selectedGenes){expressionVariance = 0}
			}
			// loop through all current selected genes and get the highest position, increment by 1 and assign position to the latest gene
			highestPos = -1
			for (key in selectedGenes){highestPos = Math.max(highestPos, selectedGenes[key])}
			highestPos += 1
			selectedGenes[gene_name] = highestPos;
			//updateGeneButtons();
			if (partyStarted){
				allocateRowAndColumnsPositions()	
			}
		}
		
		// this function is called whenever a change to the list of selected genes was made, rather than 
		// having a separate function for inserting new buttons and a function for removing buttons
		function updateGeneButtons(){
			geneMenu.innerHTML = ''
			var divContainers = []
			for (key in selectedGenes){
				var genePanel    = document.createElement('span')
					geneSelector = document.createElement('input'),
					geneRemoval  = document.createElement('input');
				geneSelector.type = 'text'
				geneSelector.value = key;
				geneSelector.gene_tag = key;
				geneSelector.onkeypress = function(event){
					if (event.which == 13){
						var enteredValue = event.target.value;
						if (gene_names.indexOf(enteredValue) != - 1){
							if (!(enteredValue in selectedGenes)){ // case gene not in selectedGens already
								// update the selectedGenes array
								var newSelectedGenes = []
								for(key in selectedGenes){
									if (key != event.target.gene_tag){
										newSelectedGenes[key] = selectedGenes[key]
									}else{
										newSelectedGenes[enteredValue] = selectedGenes[key]
									}
								}
								// update the gene buttons and the plot
								selectedGenes = newSelectedGenes
								//updateGeneButtons()
								draw()
							}else{ // case gene already in selected genes
								// generate a message for the user
								
								// restore value
								event.target.value = event.target.gene_tag
							}
						}else{
							// gene names was not found. a message will be created for the user to know
							
							// restore value
							event.target.value = event.target.gene_tag
						}
					}
				}
				geneRemoval.type = 'button'
				geneRemoval.value = 'Remove'
				geneRemoval.gene_tag = key
				geneRemoval.onclick = function(event){
					removeGeneButton(event.target.gene_tag)
				}
				genePanel.appendChild(geneSelector)
				genePanel.appendChild(geneRemoval);
				divContainers.push(genePanel)
			}
			for (key in selectedGenes){geneMenu.appendChild(divContainers[selectedGenes[key]])}
		}
		
		function removeGeneButton(gene_to_be_removed){
			var disjointPos = selectedGenes[gene_to_be_removed],
				newSelectedGenes = [];
			for (key in selectedGenes){
				if(key != gene_to_be_removed){
					newSelectedGenes[key] = selectedGenes[key]
				}
			}
			for (key in newSelectedGenes){
				if(newSelectedGenes[key] > disjointPos){
					newSelectedGenes[key] -= 1
				}
			}
			selectedGenes = newSelectedGenes;
			//updateGeneButtons()
			allocateRowAndColumnsPositions()	
		}
	
		function draw(){
			// update the output
			output_genes = ''
			for (selectedGene in selectedGenes){
				output_genes = output_genes + selectedGene + ", "
			}
			output_genes = output_genes.substring(0, output_genes.length - 2);
			
			output_cells = ''
			for (selectedCellType in selectedCellTypes){
				output_cells = output_cells + selectedCellType + ', '
			}
			output_cells = output_cells.substring(0, output_cells.length - 2);
			
			output_div_cells.innerHTML = output_cells
			output_div_genes.innerHTML = output_genes
			
			// first clear the canvas
			context.fillStyle = '#ffffff' // should change to white at the end
			context.fillRect(0, 0, canvas.width, canvas.height)
			// set parameters for label writting
			context.textAlign = 'right'
			context.textBaseline = 'middle'
			context.fillStyle = 'black'
			context.font = labelsFont;
			context.save()
			// draw if layout is cells on x axis
			if (layout == 'xcells'){
				// write the xlabels (i.e. cell names)
				for(key in selectedCellTypes){
					var col_position = selectedCellTypes[key];
					col_position = col_coordinates[col_position];
					context.translate(col_position + tileDim / 2, canvasH - padY + 6)
					context.rotate( - .23 * Math.PI)
					context.fillText(key, 0, 0)
					context.resetTransform()
				}
				// write the ylabels (i.e. gene names)
				for (key in selectedGenes){
					var row_position = selectedGenes[key]
					row_position = row_coordinates[row_position]
					context.translate(padX - 3, row_position + tileDim / 2)
					context.fillText(key, 0, -5)
					context.resetTransform()
				}
			}
			// draw if layout is cells on y axis
			else if(layout = 'xgenes'){
				for (key in selectedCellTypes){
					var row_position = selectedCellTypes[key]
					row_position = row_coordinates[row_position]
					context.translate(padX - 3, row_position + tileDim / 2)
					context.fillText(key, 0, -5)
					context.resetTransform()
				}
				for (key in selectedGenes){
					var col_position = selectedGenes[key];
					col_position = col_coordinates[col_position];
					context.translate(col_position + tileDim / 2, canvasH - padY + 6)
					context.rotate( - .23 * Math.PI)
					context.fillText(key, 0, 0)
					context.resetTransform()
				}
			}
			// draw the data viz part
			for (var cell_name in selectedCellTypes){
				var cell_name_pos = selectedCellTypes[cell_name]
				for (var gene_name in selectedGenes){
					var gene_name_pos = selectedGenes[gene_name]
					var expressionVal = getExpressionValue(cell_name, gene_name)
					if (layout == 'xcells'){
						var x_coord = selectedCellTypes[cell_name],
							y_coord = selectedGenes[gene_name];
					}else{
						var y_coord = selectedCellTypes[cell_name],
							x_coord = selectedGenes[gene_name];
					}
					x_coord = col_coordinates[x_coord] + tileDim / 2
					y_coord = row_coordinates[y_coord] + tileDim / 2
					if (plotType == 'spotplot'){
						context.fillStyle = colorCodeValue(expressionVal)
						context.beginPath()
						context.arc(x_coord, y_coord, 5 * Math.sqrt(expressionVal), 0, 2 * Math.PI, false)
						context.fill()
						context.closePath()
					}else{
						context.strokeStyle = 'black'
						context.fillStyle = colorCodeValue(expressionVal)
						context.lineWidth = 0.5;
						context.beginPath()
						context.rect(x_coord - tileDim / 2 + .5, y_coord - tileDim / 2 + .5, tileDim, tileDim)
						context.closePath()
						context.fill()
						context.stroke()
					}
				}
			}
		}
		
		function getExpressionValue(cell_name, gene_name){
			cell_name_index = cell_names.indexOf(cell_name)
			return expression_data[gene_name][cell_name_index]
		}
		
		function colorCodeValue(val){
			if (val < 4.5){
				var r = colorStart.r * (4.5 - val)/4.5 + colorMiddle.r * (val) / 4.5,
					g = colorStart.g * (4.5 - val)/4.5 + colorMiddle.g * (val) / 4.5,
					b = colorStart.b * (4.5 - val)/4.5 + colorMiddle.b * (val) / 4.5;
			}else{
				var r = colorMiddle.r * (9 - val) / 4.5 + colorStop.r * (val - 4.5) / 4.5,
					g = colorMiddle.g * (9 - val) / 4.5 + colorStop.g * (val - 4.5) / 4.5,
					b = colorMiddle.b * (9 - val) / 4.5 + colorStop.b * (val - 4.5) / 4.5;
			}
			r = parseInt(r); g = parseInt(g); b = parseInt(b);
			return 'rgb(' + r + ', ' + g + ', '+ b + ')';
		}
		
		// function that takes an event as input and return x, y values of mouse cursor
  		function getEventCoordinates(event){
  			var canvasRect = canvas.getBoundingClientRect(),
  				X = event.clientX - canvasRect.x,
  				Y = event.clientY - canvasRect.y;
  			return [X, Y]
  		}
  		
  		// function call for row dragging by mouse movement
  		function dragRow(event){
  			var coordinates = getEventCoordinates(event);
  			mouseYCurrent = coordinates[1]
  			row_coordinates[choosenRowIndex] = selectedRowPos - mouseYStart + mouseYCurrent;
  			draw()
  			
  		}
  		
  		function stopRowDragging(event){
  			canvas.removeEventListener('mousemove', dragRow)
  			canvas.removeEventListener('mouseup', stopRowDragging)
  			dragRow(event)
  			if (layout == 'xcells'){
  				selectedGenes = sortArraybyAnotherArray(selectedGenes, row_coordinates)
  			}else{
  				selectedCellTypes = sortArraybyAnotherArray(selectedCellTypes, row_coordinates)
  			}
  			allocateRowAndColumnsPositions()
  			//updateGeneButtons()
  		}
  		
  		function dragCol(event){
  			var coordinates = getEventCoordinates(event);
  			mouseXCurrent = coordinates[0]
  			col_coordinates[choosenColIndex] = selectedColPos - mouseXStart + mouseXCurrent;
  			draw()
  		}
  		
  		function stopColDragging(event){
  			canvas.removeEventListener('mousemove', dragCol)
  			canvas.removeEventListener('mouseup', stopColDragging)
  			if (layout == 'xgenes'){
  				selectedGenes = sortArraybyAnotherArray(selectedGenes, col_coordinates)
  			}else{
  				selectedCellTypes = sortArraybyAnotherArray(selectedCellTypes, col_coordinates)
  			}
  			allocateRowAndColumnsPositions()
  			//updateGeneButtons()
  		}
  		
  		// sorts a named array by a second array of values
  		function sortArraybyAnotherArray(namedArr, valArr){
  			var names        = [],
  				newNamedArr  = [],
  				sortedValues = valArr.slice()
  			sortedValues.sort(function(a, b){return a - b})
  			for(key in namedArr){names.push(key)}
  			for (var i=0; i<valArr.length; i++){
  				var nameIndex = valArr.indexOf(sortedValues[i])
  				newNamedArr[names[nameIndex]] = i
  			}
  			return newNamedArr
  		}
		
		// instructions to be execute upon complete page load
		window.onload = function(){
			// fire procedure for making the cell menu
			cell_names.forEach(function(val,i){
				var cell_span  = document.createElement('span'),
					cell_input = document.createElement('input'),
					tag_span   = document.createElement('span');
				tag_span.innerHTML = val
				tag_span.style.paddingRight = '.8em'
				cell_input.type = 'checkbox'
				cell_input.checked = true
				cell_input.value = val
				cell_input.onchange = updateSelectedCellTypes
				cell_span.appendChild(cell_input)
				cell_span.appendChild(tag_span)
				cellMenu.appendChild(cell_span)
			})
			
			// initialise selectedGenes and selectedCellTypes from starting values:
			let initialCellTypeNames = Array.from(document.querySelectorAll("cell-menu input[type='checkbox']"), c => c.value);
			for (var i=0;i<initialCellTypeNames.length;i++){
				selectedCellTypes[initialCellTypeNames[i]] = i
			}

			let initialGeneSymbols = verifyData(document.querySelector("gene-selector").items);
			for (var i=0;i<initialGeneSymbols.length;i++){
				selectedGenes[initialGeneSymbols[i]] = i
			}

			// everything should be ready to get the interactivity started
			partyStarted = true
			allocateRowAndColumnsPositions()


			// add event listeners
			canvas.addEventListener('mousedown', function(event){
				var coordinates = getEventCoordinates(event),
					clickX      = coordinates[0],
					clickY      = coordinates[1];
				// if selection is in rows
				if (clickX < padX && clickY < (canvasH - padY)){
					var rowDists = row_coordinates.map(function(val){return val + tileDim / 2})
					               .map(function(val){return Math.abs(clickY - val)}),
						smallestDist = rowDists.reduce(function(a, b){return Math.min(a, b)});
					// store chosen row index and Y position of click event firing
					choosenRowIndex = rowDists.indexOf(smallestDist);
					mouseYStart     = clickY;
					selectedRowPos  = row_coordinates[choosenRowIndex]
					// add event listener for mouse movement
					canvas.addEventListener('mousemove', dragRow)
					canvas.addEventListener('mouseup', stopRowDragging);
				}else if(clickY > (canvasH - padY)){ // check is selection is in columns
					var trueX   = clickX + (clickY - canvasH + padY) * Math.tan(Math.PI / 2 -.23 * Math.PI),
						colDist = col_coordinates.map(function(val){return val + tileDim / 2})
												 .map(function(val){return Math.abs(trueX - val)})
						smallestDist = colDist.reduce(function(a, b){return Math.min(a, b)})
					choosenColIndex = colDist.indexOf(smallestDist)
					withinLabel = true; // if click is within label at choosen index
					if(withinLabel){
						mouseXStart = clickX
						selectedColPos = col_coordinates[choosenColIndex]
						canvas.addEventListener('mousemove', dragCol)
						canvas.addEventListener('mouseup', stopColDragging);
					}
				}
			})

			// add event listener to cell type gene list reader 
			cellTypeListReader.onkeypress = function(event){
				if (event.which == 13){
					var cellType_list_input = cellTypeListReader.value
					cellTypeListReader.value = '';
					cellType_list_input = cellType_list_input.replace(/, /g, 'splitpoint')
					                                         .replace(/,/g, 'splitpoint')
					                                         .replace(/; /g, 'splitpoint')
					                                         .replace(/;/g, 'splitpoint')
					                                         .replace(/\t /g, 'splitpoint')
					                                         .replace(/\t/g, 'splitpoint')
					                                         .split('splitpoint');
					// remove all the empty strings in cell types lists
					new_cellType_list_input = []
					for(var i=0;i<cellType_list_input.length;i++){
						if(cellType_list_input[i] != ''){
							new_cellType_list_input.push(cellType_list_input[i])
						}
					}
					
					cellType_list_input = new_cellType_list_input;
					// check if all gene in the list exist in the data
					var allInList = true,
						cellType_not_found = [];
					for(var i=0;i<cellType_list_input.length;i++){
						if (cell_names.indexOf(cellType_list_input[i]) == -1){
							allInList = false;
							cellType_not_found.push(cellType_list_input[i])
						}
					}
					// check for repeats
					var counter_array = []
					for (var i=0;i<cellType_list_input.length;i++){
						celltype_name = cellType_list_input[i]
						if (celltype_name in counter_array){
							counter_array[celltype_name] += 1
						}else {
							counter_array[celltype_name] = 1
						}
					}
					var duplicates      = [],
						duplicatesExist = false;
					for(key in counter_array){if(counter_array[key] > 1){duplicates.push(key); duplicatesExist = true}}
					if (duplicatesExist){
						msg = ''
						for (var i=0; i < duplicates.length; i++){
							msg = msg + duplicates[i] + ', '
						}
						msg = 'There are duplicates: ' + msg
						cellTypeListReader.value = msg;
						return 1;
					}
					// if all are in list and the are not repeats render the plot and make the gene buttons
					if (allInList && cellType_list_input.length > 0){
						selectedCellTypes = []
						for (var i=0;i<cellType_list_input.length;i++){
							selectedCellTypes[cellType_list_input[i]] = i
						}
						allocateRowAndColumnsPositions()
						//updateGeneButtons()
						// update the cell types radio buttons
						cell_buttons = cellMenu.children;
						for(var i=0;i<cell_buttons.length;i++){
							var spanID = cell_buttons[i].children[0];
							if(cellType_list_input.indexOf(cellMenu.children[i].children[0].value) >= 0 ){
								spanID.checked = true;
							}else{
								spanID.checked = false;
							}
						}
					}else{
						var msg = '';
						for (i=0;i<cellType_not_found.length;i++){
							msg = msg + cellType_not_found[i] + ', '
						}
						msg = 'Cell types not found: ' + msg;
						cellTypeListReader.value = msg;
					}
					// if nothing was passed:
					if(cellType_list_input.length == 0){
						cellTypeListReader.value = ''
					}
				}
			}

			// add event listener to the gene list reader text area so that it can take list of genes as input for selected genes
			document.querySelector("gene-selector").addEventListener("update", function(u) {
				var gene_list_input  = this.items.slice();

				gene_list_input = gene_list_input.filter(gene => Boolean(expression_data[gene])); // exclude values not in the dataset
				gene_list_input = gene_list_input.filter(gene => Boolean(gene)); // exclude missing & blank values (i.e.
				gene_list_input = gene_list_input.filter((gene, i, list) => list.indexOf(gene) === i) // remove duplicate values N.B. the gene-selector component doesn't allow duplicates, but this logic could be reused with other components that don't enforce that.

				// by this point we have a cleaned list of unique genes that are in the dataset to display on the chart.

				// if there are existing genes we don't want to re-order the chart, so we order and filter appropriately:
				if (Object.keys(selectedGenes).length > 0) {
					// get a sorted list of existing gene names in selectedGenes;
					let currGenes = Object.keys(selectedGenes.sort((a,b) => a-b));
					// filter only genes that are in the input list (i.e. remove any genes that aren't)
					currGenes = currGenes.filter(gene => gene_list_input.includes(gene))
					// concatenate with genes from input that *aren't* in selected genes:
					currGenes = currGenes.concat(gene_list_input.filter(gene => !currGenes.includes(gene)));

					gene_list_input = currGenes;
				}

				// clear selected genes then reallocate new positions from input:
				selectedGenes = [];
				for (var i=0;i<gene_list_input.length;i++){
					selectedGenes[gene_list_input[i]] = i
				}
				
				// update chart
				allocateRowAndColumnsPositions()
			});
		}
	</script>
	<div style="float: clear;""><hr><span style="font-size:0.8em;">This data portal was created using the interactive_heatmap_dotplot tool (<a href="https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data">github link</a>) developed by Dorin-Mirel Popescu and Issac Goh, maintained by Issac Goh</span><hr></div>
</body>
</html>
