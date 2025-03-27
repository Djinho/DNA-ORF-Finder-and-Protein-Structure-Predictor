Sys.setenv(RETICULATE_PYTHON = "C:/Users/Djinh/conda/envs/maker_env_windows/python.exe")
library(shiny)
library(bslib)
library(reticulate)

# Define codon-to-amino acid map based on the genetic code table
codon_table <- list(
  UUU = "Phe", UUC = "Phe", UUA = "Leu", UUG = "Leu",
  CUU = "Leu", CUC = "Leu", CUA = "Leu", CUG = "Leu",
  AUU = "Ile", AUC = "Ile", AUA = "Ile", AUG = "Met",
  GUU = "Val", GUC = "Val", GUA = "Val", GUG = "Val",
  UCU = "Ser", UCC = "Ser", UCA = "Ser", UCG = "Ser",
  CCU = "Pro", CCC = "Pro", CCA = "Pro", CCG = "Pro",
  ACU = "Thr", ACC = "Thr", ACA = "Thr", ACG = "Thr",
  GCU = "Ala", GCC = "Ala", GCA = "Ala", GCG = "Ala",
  UAU = "Tyr", UAC = "Tyr", UAA = "Stop", UAG = "Stop",
  CAU = "His", CAC = "His", CAA = "Gln", CAG = "Gln",
  AAU = "Asn", AAC = "Asn", AAA = "Lys", AAG = "Lys",
  GAU = "Asp", GAC = "Asp", GAA = "Glu", GAG = "Glu",
  UGU = "Cys", UGC = "Cys", UGA = "Stop", UGG = "Trp",
  CGU = "Arg", CGC = "Arg", CGA = "Arg", CGG = "Arg",
  AGU = "Ser", AGC = "Ser", AGA = "Arg", AGG = "Arg",
  GGU = "Gly", GGC = "Gly", GGA = "Gly", GGG = "Gly"
)

# Function to translate mRNA sequence to amino acids
translate_to_amino_acids <- function(mrna) {
  amino_acids <- ""
  for (i in seq(1, nchar(mrna) - 2, by = 3)) {
    codon <- substr(mrna, i, i + 2)
    if (codon_table[[codon]] == "Stop") break
    amino_acids <- paste0(amino_acids, codon_table[[codon]], "-")
  }
  return(substr(amino_acids, 1, nchar(amino_acids) - 1)) # Remove trailing "-"
}

# Function to find ORFs in a DNA sequence
find_orfs <- function(sequence) {
  start_codon <- "ATG"
  stop_codons <- c("TAA", "TAG", "TGA")
  orfs <- list()
  
  sequence <- toupper(gsub("[^ATCG]", "", sequence))  # Clean and validate DNA sequence
  
  for (frame in 0:2) {
    orf <- ""
    started <- FALSE
    for (i in seq(frame + 1, nchar(sequence) - 2, by = 3)) {  # Adjusted to nchar(sequence) - 2
      codon <- substr(sequence, i, i + 2)
      
      if (!started && codon == start_codon) {
        orf <- codon
        started <- TRUE
      } else if (started) {
        orf <- paste0(orf, codon)
        if (codon %in% stop_codons) {
          orfs <- append(orfs, orf)
          orf <- ""
          started <- FALSE
        }
      }
    }
  }
  
  if (length(orfs) == 0) {
    return("No open reading frames (ORFs) found.")
  } else {
    return(orfs)
  }
}

# Define the UI with the Minty theme
ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
  
  titlePanel("DNA ORF Finder and Protein Structure Predictor"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput("dnaSeq", "Enter DNA Sequence", 
                    placeholder = "Enter a DNA sequence in ATCG format..."),
      actionButton("findOrfs", "Find ORFs and Translate"),
      hr(),
      textAreaInput("proteinSeq", "Enter Protein Sequence for Secondary Structure Prediction",
                    placeholder = "Enter a protein sequence in FASTA format..."),
      actionButton("predictStructure", "Predict Secondary Structure")
    ),
    
    mainPanel(
      h3("Open Reading Frames (ORFs) Found"),
      verbatimTextOutput("orfOutput"),
      hr(),
      h3("Predicted Secondary Structure"),
      uiOutput("structureOutput")
    )
  )
)

# Define server logic for ORF finding, translation, and structure prediction
server <- function(input, output, session) {
  
  # Display ORFs when the button is clicked
  output$orfOutput <- renderPrint({
    req(input$findOrfs)  # Ensure button is clicked
    dna_seq <- input$dnaSeq
    if (nchar(dna_seq) == 0) {
      return("Please enter a DNA sequence.")
    }
    
    # Find ORFs in the DNA sequence
    orfs <- find_orfs(dna_seq)
    if (is.character(orfs) && orfs == "No open reading frames (ORFs) found.") {
      return(orfs)
    }
    
    # Translate each ORF to an amino acid sequence
    translated_orfs <- list()
    for (orf in orfs) {
      mrna <- gsub("T", "U", orf)  # Convert DNA to mRNA
      amino_acids <- translate_to_amino_acids(mrna)
      translated_orfs[[orf]] <- amino_acids
    }
    
    return(translated_orfs)
  })
  
  # Predict secondary structure when the button is clicked
  output$structureOutput <- renderUI({
    req(input$predictStructure)  # Ensure button is clicked
    protein_seq <- input$proteinSeq
    if (nchar(protein_seq) == 0) {
      return(HTML("<p>Please enter a protein sequence.</p>"))
    }
    
    # Write the protein sequence to a temporary FASTA file
    temp_fasta <- tempfile(fileext = ".fas")
    writeLines(protein_seq, temp_fasta)
    
    # Run S4PRED on the temporary FASTA file
    result <- system(paste(Sys.getenv("RETICULATE_PYTHON"), "C:/Users/Djinh/s4pred/run_model.py", temp_fasta), intern = TRUE)
    
    # Filter out unnecessary warnings and display predictions clearly
    clean_result <- result[!grepl("FutureWarning", result)]  # Remove warnings
    clean_result <- gsub("^\\s+|\\s+$", "", clean_result)  # Trim whitespace
    
    # Only display structured prediction output
    prediction_start <- grep("PSIPRED VFORMAT", clean_result)
    if (length(prediction_start) > 0) {
      formatted_output <- ""
      sequence_counter <- 1
      
      for (i in prediction_start) {
        if (sequence_counter > 1) {
          # Add separator for each new sequence
          formatted_output <- paste0(formatted_output, "<hr style='border-top: 2px solid #ccc;'>")
        }
        
        sequence_data <- clean_result[i:length(clean_result)]
        formatted_rows <- lapply(sequence_data, function(line) {
          columns <- unlist(strsplit(line, "\\s+"))
          if (length(columns) == 6) {
            return(
              paste0("<tr><td>", columns[1], "</td><td>", columns[2], "</td><td>", columns[3], "</td><td>", 
                     columns[4], "</td><td>", columns[5], "</td><td>", columns[6], "</td></tr>")
            )
          }
          return(NULL)
        })
        
        table_html <- paste(
          "<table style='width:100%; border-collapse:collapse;'>",
          "<tr><th>Position</th><th>Amino Acid</th><th>Structure</th><th>Coil Conf.</th><th>Helix Conf.</th><th>Strand Conf.</th></tr>",
          paste(unlist(formatted_rows), collapse = ""),
          "</table>"
        )
        
        formatted_output <- paste(formatted_output, table_html)
        sequence_counter <- sequence_counter + 1
      }
      
      return(HTML(formatted_output))
    } else {
      return(HTML("<p>No prediction available.</p>"))
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)

