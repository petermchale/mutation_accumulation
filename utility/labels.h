#ifndef LABELS_H
#define	LABELS_H

/*************************************************************************/

namespace labels {

    /** 
     * create label appropriate for mean 
     */
    template <class sample_type>
    const std::string create_label(probability::Mean<sample_type>) {

        return std::string("mean");
    }


}

#endif	/* LABELS_H */

