#ifndef STRINGS_H
#define	STRINGS_H

#include <string> 
#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp> // boost::split
#include <boost/lexical_cast.hpp> //boost::lexical_cast

#include <mutation_accumulation/array/array2D.h> // array::Array2D 

/*************************************************************************/

/**
 * useful functions for manipulating strings\n
 *
 * Because templates are compiled when required, \n
 * the implementation of a template function (or class) \n
 * must be in the same file as its declaration
 */
namespace strings {

    /** 
     * read input file stream into a vector of strings 
     */
    const std::vector<std::string> read_vectorStrings(std::ifstream &ifs, const std::string &delimiter = ",") {

        std::string text_line;
        std::getline(ifs, text_line);

        std::vector<std::string> tokens;
        boost::split(tokens, text_line, boost::is_any_of("="));

        const std::string values = tokens.at(1);
        boost::split(tokens, values, boost::is_any_of(delimiter));

        for (int ii = 0; ii < tokens.size(); ii++) 
            boost::trim(tokens.at(ii));

        return tokens;

    }

    /** 
     * convert vector of strings to vector of values 
     */
    template <class parameter_type>
    const std::vector<parameter_type> vectorStrings_to_vectorValues(const std::vector<std::string> &parameter_strings) {

        std::vector<parameter_type> parameter_values(parameter_strings.size(), static_cast<parameter_type> (0));

        for (int ii = 0; ii < parameter_strings.size(); ii++) {

            // cast first to <double> in case token is written in exponential notation
            const double parameter_value_double = boost::lexical_cast<double> (boost::trim_copy(parameter_strings.at(ii)));
            
            // then cast to parameter_type
            parameter_values.at(ii) = static_cast<parameter_type> (parameter_value_double);
        }

        return parameter_values;

    }

    namespace strings_detail {

        template <class parameter_type>
        const parameter_type do_parse_scalar(std::ifstream &ifs) {

            std::string text_line;
            std::getline(ifs, text_line);

            std::vector<std::string> tokens;

            boost::split(tokens, text_line, boost::is_any_of("="));

            // cast first to <double> in case token is written in exponential notation
            const double parameter_value_double = boost::lexical_cast<double> (boost::trim_copy(tokens.at(1)));

            // then cast to parameter_type
            return static_cast<parameter_type> (parameter_value_double);

        }

        template <class parameter_type>
        const std::vector<parameter_type> do_parse_vector(std::ifstream &ifs, const std::string &delimiter = ",") {

            return vectorStrings_to_vectorValues<parameter_type>(read_vectorStrings(ifs, delimiter));
        }

        template <class parameter_type>
        const array::Array2D<parameter_type> do_parse_matrix(std::ifstream &ifs) {

            std::string line_string;
            std::getline(ifs, line_string);

            std::vector<std::string> line_tokens;
            boost::split(line_tokens, line_string, boost::is_any_of("="));

            const std::string matrix_string = line_tokens.at(1);
            std::vector<std::string> row_tokens;
            boost::split(row_tokens, matrix_string, boost::is_any_of(";"));

            std::vector<std::vector<parameter_type> > matrix(row_tokens.size());

            for (int row = 0; row < row_tokens.size(); row++) {

                const std::string row_string = boost::trim_copy(row_tokens.at(row));
                std::vector<std::string> col_tokens;
                boost::split(col_tokens, row_string, boost::is_any_of(","));

                std::vector<parameter_type> vector(col_tokens.size());

                for (int col = 0; col < col_tokens.size(); col++) {

                    // cast first to <double> in case token is written in exponential notation
                    const double parameter_value_double = boost::lexical_cast<double> (boost::trim_copy(col_tokens.at(col)));

                    // then cast to parameter_type
                    vector.at(col) = static_cast<parameter_type> (parameter_value_double);
                }

                matrix.at(row) = vector;

            }

            /* check that number of columns supplied by user is identical for each row */
            for (int row = 0; row < matrix.size(); row++)
                assert(matrix.at(row).size() == matrix.at(0).size());

            /* build Array2D object */
            const int number_rows = matrix.size();
            const int number_cols = matrix.at(0).size();
            array::Array2D<parameter_type> array2D(number_rows, number_cols);
            for (int row = 0; row < number_rows; row++)
                for (int col = 0; col < number_cols; col++)
                    array2D.at(row, col) = matrix.at(row).at(col);

            /* return Array2D object */
            return array2D;

        }

    }

    /**
     * find value of parameter in next line of text\n
     * text takes the form "parameter = value"\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    void parse_scalar(std::ifstream &ifs, parameter_type & parameter_value) {

        parameter_value = strings_detail::do_parse_scalar<parameter_type > (ifs);

    }

    /**
     * find value of parameter in next line of text\n
     * text takes the form "parameter = value"\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    const parameter_type parse_scalar(std::ifstream & ifs) {

        return strings_detail::do_parse_scalar<parameter_type > (ifs);
    }

    /**
     * find values of parameter in next line of text\n
     * text takes the form "parameter = value1, value2, .. ,valueN"\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    void parse_vector(std::ifstream &ifs, std::vector<parameter_type> &parameter_values, const std::string &delimiter = ",") {

        parameter_values = strings_detail::do_parse_vector<parameter_type > (ifs, delimiter);

    }

    /**
     * find values of parameter in next line of text\n
     * text takes the form "parameter = value1, value2, .. ,valueN"\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    const std::vector<parameter_type> parse_vector(std::ifstream &ifs, const std::string &delimiter = ",") {

        return strings_detail::do_parse_vector<parameter_type > (ifs, delimiter);

    }

    /**
     * find values of parameter in next line of text\n
     * text takes the form "parameter = value1a, value1b, ...; value2a, value2b, ... ;valueNa, valueNb, ..."\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    void parse_matrix(std::ifstream &ifs, array::Array2D<parameter_type> &parameter_values) {

        parameter_values = strings_detail::do_parse_matrix<parameter_type > (ifs);

    }

    /**
     * find values of parameter in next line of text\n
     * text takes the form "parameter = value1a, value1b, ...; value2a, value2b, ... ;valueNa, valueNb, ..."\n
     * parameter_type can be built-in or user-defined
     */
    template <class parameter_type>
    const array::Array2D<parameter_type> parse_matrix(std::ifstream & ifs) {

        return strings_detail::do_parse_matrix<parameter_type > (ifs);

    }



}


#endif	/* STRINGS_H */

