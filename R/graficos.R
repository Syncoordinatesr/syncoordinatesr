
#' @title  An auxiliary function to generate traceplot graphs based on MCMC results for further analysis.
#'
#' @name  graphs_traceplots
#'
#' @param resultados A list containing MCMC results, including elements like "S", "burn", "alfa", "mu", "tau.theta", "tau.phi", and "tau.e".
#' @param dados A list containing dataset information.
#' @param plot_color A character indicating the color for the plot lines.
#' @param add_trendline A logical indicating whether to add a trendline to the plots.
#' @param save_file A logical indicating whether to save the generated plots to files.
#'
#' @return It generates traceplot graphs based on the provided MCMC results and saves them if specified.
#'
#' @examples
#' # Example usage:
#' # graphs_traceplots(resultados, dados, plot_color = "blue", add_trendline = TRUE, save_file = TRUE)
#'
#' @references
#' NUNES, Letícia. Métodos de Simulação de Dados Geográficos Sintéticos Para Bases Confidenciais. *Dissertação de Mestrado*, [s. l.], 2018.
#' Disponível em: \url:{http//est.ufmg.br/portal/arquivos/mestrado/dissertacoes/dissertacao_Leticia_Silva_Nunes.pdf}. Acesso em: 2 mar. 2022.
#'
#' @export
#'


# Traceplots Function -----------------------------------------------------
graphs_traceplots  <- function(resultados, dados, plot_color, add_trendline, save_file) {

      num <- length(resultados)
      if (
        !is.list(resultados) ||  # Ensures resultados as a list
        !is.list(dados) ||  # Ensures dados as a list
        num < 4||
        !is.logical(save_file) || # Ensures save_file is a logical argument
        !is.character(plot_color)||   # Ensures plot_color is a character
        !is.logical(add_trendline)  # Ensures add_trendline is a logical argument

      ) {
        stop("Arguments are not in the correct format. Please ensure that results and data are provided as lists,
          save file and add_trendline are specified as T or F, and check that plot_color is a character.")
      }

      else{

        # Extracting specific elements from the 'resultados' list
        S = resultados[["S"]]
        burn = resultados[["burn"]]
        alfa = resultados[["alfa"]]
        mu = resultados[["mu"]]
        tau.theta = resultados[["tau.theta"]]
        tau.phi = resultados[["tau.phi"]]
        tau.e = resultados[["tau.e"]]

        # Calculating the number of rows in 'alfa' and 'tau.phi' matrices
        num_alfa <- nrow(alfa)
        num_tau_phi <- nrow(tau.phi)

        # Calculating the total number of combinations (rows in 'alfa' and 'tau.phi' matrices plus 3)
        num_combinacoes <- num_alfa + num_tau_phi + 3

        # Calculating the number of rows and columns for the plot layout
        nrows <- ceiling(sqrt(num_combinacoes))
        ncols <- ceiling(num_combinacoes / nrows)

        par(mfrow = c(nrows, ncols),oma = c(1,1,1,1),mar = c(1,1,2,1))

        if (add_trendline == T) { # Plot graphs with central tendency line

          # Mu
          plot(mu[1:S], type = "l", xlab = "", ylab = "", main= "(Intercept)",col = plot_color,yaxt = "n",xaxt = "n")
          abline(lm(mu[1:S] ~ seq_along(mu[1:S])), col = "black",lty = 2)

          # Tau.theta
          plot(tau.theta[1:S], type = "l", xlab = "", ylab = "", main = "Tau.theta",col = plot_color,yaxt = "n",xaxt = "n")
          abline(lm(tau.theta[1:S] ~ seq_along(tau.theta[1:S])), col = "black",lty = 2)

          # Tau.e
          plot(tau.e[1:S], type = "l", xlab = "", ylab = "", main = "Tau.e",col = plot_color,yaxt = "n",xaxt = "n")
          abline(lm(tau.e[1:S] ~ seq_along(tau.e[1:S])), col = "black",lty = 2)


          # Alfa
          for (i in 1:num_alfa) {
            plot(alfa[i, 1:S], type = "l", xlab = "", ylab = "", main= paste("Alfa.", i),col = plot_color,yaxt = "n",xaxt = "n")
            abline(lm(alfa[i, 1:S] ~ seq_along(alfa[i, 1:S])), col = "black",lty = 2)


          }

          # Tau.phi
          for (i in 1:num_tau_phi) {
            plot(tau.phi[i, 1:S], type = "l", xlab = "", ylab = "", main= paste("Tau.phi.", i),col = plot_color,yaxt = "n",xaxt = "n")
            abline(lm(tau.phi[i, 1:S] ~ seq_along(tau.phi[i, 1:S])), col = "black",lty = 2)

          }

          par(mfrow = c(1, 1))

        }

      else{

        # Mu
        plot(mu[1:S], type = "l", xlab = "", ylab = "", main= "(Intercept)",col = plot_color,yaxt = "n",xaxt = "n")

        # Tau.theta
        plot(tau.theta[1:S], type = "l", xlab = "", ylab = "", main = "Tau.theta",col = plot_color,yaxt = "n",xaxt = "n")

        # Tau.e
        plot(tau.e[1:S], type = "l", xlab = "", ylab = "", main = "Tau.e",col = plot_color,yaxt = "n",xaxt = "n")

        # Alfa
        for (i in 1:num_alfa) {
          plot(alfa[i, 1:S], type = "l", xlab = "", ylab = "", main= paste("Alfa.", i),col = plot_color,yaxt = "n",xaxt = "n")

        }

        # Tau.phi
        for (i in 1:num_tau_phi) {
          plot(tau.phi[i, 1:S], type = "l", xlab = "", ylab = "", main= paste("Tau.phi.", i),col = plot_color,yaxt = "n",xaxt = "n")
        }

        par(mfrow = c(1, 1))
      }

      }


      if (save_file == T) {

        # Creating a folder to store the traceplot graphs
        graphs_traceplots <- "graphs/traceplots"
        if (!file.exists(graphs_traceplots)) {
          dir.create(graphs_traceplots, recursive = TRUE)
        }

        if (add_trendline == T) {

          #Mu
          pdf_mu <- file.path(graphs_traceplots, "plot_mu.pdf")
          pdf(file = pdf_mu)
          plot(mu[1:S], type = "l", xlab = "Iterations", ylab = "mu", col= plot_color)
          abline(lm(mu[1:S] ~ seq_along(mu[1:S])), col = "black",lty = 2)
          dev.off()


          #Tau.theta
          pdf_tau_theta <- file.path(graphs_traceplots, "plot_tau_theta.pdf")
          pdf(file = pdf_tau_theta)
          plot(tau.theta[1:S], type = "l", xlab = "Iterations", ylab = "tau.theta", col= plot_color)
          abline(lm(tau.theta[1:S] ~ seq_along(tau.theta[1:S])), col = "black",lty = 2)
          dev.off()

          #Tau.e
          pdf_tau_e <- file.path(graphs_traceplots, "plot_tau_e.pdf")
          pdf(file = pdf_tau_e)
          plot(tau.e[1:S], type = "l", xlab = "Iterations", ylab = "tau.e",col= plot_color)
          abline(lm(tau.e[1:S] ~ seq_along(tau.e[1:S])), col = "black",lty = 2)
          dev.off()

          #Alfa
          num_alfa <- nrow(alfa)
          for (i in 1:num_alfa) {
            pdf_alfa <- file.path(graphs_traceplots, paste0("plot_alfa", i, ".pdf"))
            pdf(file = pdf_alfa)
            plot(alfa[i, 1:S], type = "l", xlab = "Iterations", ylab = paste("alfa (", i, ")"), col= plot_color)
            abline(lm(alfa[i, 1:S] ~ seq_along(alfa[i, 1:S])), col = "black",lty = 2)
            dev.off()
        }

          #Tau.phi
          num_tau_phi <- nrow(tau.phi)
          for (i in 1:num_tau_phi) {
            pdf_tau_phi <- file.path(graphs_traceplots, paste0("plot_tau_phi", i, ".pdf"))
            pdf(file = pdf_tau_phi)
            plot(tau.phi[i,1:S], type = "l", xlab = "Iterations", ylab = paste("tau.phi (", i, ")"), col= plot_color)
            abline(lm(tau.phi[i, 1:S] ~ seq_along(tau.phi[i, 1:S])), col = "black",lty = 2)
            dev.off()

        }
        }
          else{

          #Mu
          pdf_mu <- file.path(graphs_traceplots, "plot_mu.pdf")
          pdf(file = pdf_mu)
          plot(mu[1:S], type = "l", xlab = "Iterations", ylab = "mu", col= plot_color)
          dev.off()


          #Tau.theta
          pdf_tau_theta <- file.path(graphs_traceplots, "plot_tau_theta.pdf")
          pdf(file = pdf_tau_theta)
          plot(tau.theta[1:S], type = "l", xlab = "Iterations", ylab = "tau.theta", col= plot_color)
          dev.off()

          #Tau.e
          pdf_tau_e <- file.path(graphs_traceplots, "plot_tau_e.pdf")
          pdf(file = pdf_tau_e)
          plot(tau.e[1:S], type = "l", xlab = "Iterations", ylab = "tau.e",col= plot_color)
          dev.off()

          #Alfa
          num_alfa <- nrow(alfa)
          for (i in 1:num_alfa) {
            pdf_alfa <- file.path(graphs_traceplots, paste0("plot_alfa", i, ".pdf"))
            pdf(file = pdf_alfa)
            plot(alfa[i, 1:S], type = "l", xlab = "Iterations", ylab = paste("alfa (", i, ")"), col= plot_color)
            dev.off()
          }

          #Tau.phi
          num_tau_phi <- nrow(tau.phi)
          for (i in 1:num_tau_phi) {
            pdf_tau_phi <- file.path(graphs_traceplots, paste0("plot_tau_phi", i, ".pdf"))
            pdf(file = pdf_tau_phi)
            plot(tau.phi[i,1:S], type = "l", xlab = "Iterations", ylab = paste("tau.phi (", i, ")"), col= plot_color)
            dev.off()

          }
        }

      }

}


